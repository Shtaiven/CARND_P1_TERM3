#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <cfloat>

using namespace std;

// for convenience
using json = nlohmann::json;

/*****************************************************************************
 * Global Constants and Structs
 *****************************************************************************/

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Holds a series of points
struct trajectory_t {
  vector<double> x_vec;
  vector<double> y_vec;
  double ref_vel;
  int final_lane;
};

struct point_t {
  double x;
  double y;
};

// Vehicle data POD
struct vehicle_t {
  double id;
  double x;
  double y;
  double vx;
  double vy;
  double s;
  double d;
  double speed;
  double yaw;
};

// Enumerate states
enum vehicle_state_t {
  KEEP_LANE,
  PREP_LANE_CHANGE_LEFT,
  LANE_CHANGE_LEFT,
  PREP_LANE_CHANGE_RIGHT,
  LANE_CHANGE_RIGHT
};

/*****************************************************************************
 * Global Function Definitions
 *****************************************************************************/

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

double distance(vehicle_t p1, vehicle_t p2)
{
  return distance(p1.x, p1.y, p2.x, p2.y);
}

void merge_trajectories(trajectory_t & t1, trajectory_t t2)
{
  t1.x_vec.reserve(t1.x_vec.size() + distance(t2.x_vec.begin(), t2.x_vec.end()));
  t1.x_vec.insert(t1.x_vec.end(), t2.x_vec.begin(), t2.x_vec.end());

  t1.y_vec.reserve(t1.y_vec.size() + distance(t2.y_vec.begin(), t2.y_vec.end()));
  t1.y_vec.insert(t1.y_vec.end(), t2.y_vec.begin(), t2.y_vec.end());
}

int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
  {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen)
    {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}

// Below are 3 functions used for line segment intersection calculations taken from https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(point_t p, point_t q, point_t r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
       return true;

    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(point_t p, point_t q, point_t r)
{
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
    // for details of below formula.
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;  // colinear

    return (val > 0)? 1: 2; // clock or counterclock wise
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(point_t p1, point_t q1, point_t p2, point_t q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

/*****************************************************************************
 * Trajectory Generation
 *****************************************************************************/

vector<vehicle_state_t> successor_states(vehicle_state_t state, int lane, int lanes_available)
{
  /*
  Provides the possible next states given the current state for the FSM
  discussed in the course, with the exception that lane changes happen
  instantaneously, so LCL and LCR can only transition back to KL.
  */
  vector<vehicle_state_t> states;
  states.push_back(KEEP_LANE);

  if (state == KEEP_LANE)
  {
    states.push_back(LANE_CHANGE_LEFT);
    states.push_back(LANE_CHANGE_RIGHT);
  }

  return states;
}

// Used for calculating trajectories of observed vehicles
trajectory_t constant_speed_trajectory(vehicle_t vehicle, int num_points,
                                       trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  trajectory_t trajectory;
  trajectory.ref_vel = vehicle.speed;
  trajectory.final_lane = (int) floor(vehicle.d/4);

  if (!vehicle.speed)
  {
    vector<double> x_vec(num_points, vehicle.x);
    vector<double> y_vec(num_points, vehicle.y);
    trajectory.x_vec = x_vec;
    trajectory.y_vec = y_vec;
    return trajectory;
  }

  vector<double> pts_x;
  vector<double> pts_y;

  // Define our starting reference as the position of the car currently
  double ref_x = vehicle.x;
  double ref_y = vehicle.y;
  double ref_yaw = vehicle.yaw;

  // Start our path by looking at the previous point traversed and the current one
  double prev_car_x = vehicle.x - cos(ref_yaw);
  double prev_car_y = vehicle.y - sin(ref_yaw);

  pts_x.push_back(prev_car_x);
  pts_x.push_back(vehicle.x);

  pts_y.push_back(prev_car_y);
  pts_y.push_back(vehicle.y);

  vector<double> next_wp0;
  vector<double> next_wp1;
  vector<double> next_wp2;

  // Add even 30 meter points in front of the car's reference position to our list of next points
  if (vehicle.speed > 0)
  {
    next_wp0 = getXY(vehicle.s+30, vehicle.d, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    next_wp1 = getXY(vehicle.s+60, vehicle.d, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    next_wp2 = getXY(vehicle.s+90, vehicle.d, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
  }
  else
  {
    next_wp0 = getXY(vehicle.s-30, vehicle.d, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    next_wp1 = getXY(vehicle.s-60, vehicle.d, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    next_wp2 = getXY(vehicle.s-90, vehicle.d, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
  }

  pts_x.push_back(next_wp0[0]);
  pts_x.push_back(next_wp1[0]);
  pts_x.push_back(next_wp2[0]);

  pts_y.push_back(next_wp0[1]);
  pts_y.push_back(next_wp1[1]);
  pts_y.push_back(next_wp2[1]);

  // Define car's position as (0,0) with yaw 0
  for (int i = 0; i < pts_x.size(); i++)
  {
    double shift_x = pts_x[i] - ref_x;
    double shift_y = pts_y[i] - ref_y;

    pts_x[i] = (shift_x * cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
    pts_y[i] = (shift_x * sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
  }

  // Create a spline curve between the points we generated
  tk::spline s;
  s.set_points(pts_x, pts_y);

  // break up spline points so that we travel at our reference velocity
  double target_x = 30.0;
  double target_y = s(target_x);
  double target_dist = sqrt(target_x*target_x+target_y*target_y);

  double x_add_on = 0;

  // Interpolate between the points we generated previously using the spline function until we have a total of 50 points
  for (int i = 0; i < num_points; i++)
  {
    double N = target_dist/(0.02*trajectory.ref_vel/2.24);
    double x_point = x_add_on+target_x/N;
    double y_point = s(x_point);

    x_add_on = x_point;

    double x_ref = x_point;
    double y_ref = y_point;

    // Transfer car coordinates back to normal instead of relative to car
    x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
    y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));

    x_point += ref_x;
    y_point += ref_y;

    trajectory.x_vec.push_back(x_point);
    trajectory.y_vec.push_back(y_point);
  }

  return trajectory;
}

trajectory_t generate_trajectory(vehicle_t vehicle, vehicle_state_t state, vector<vehicle_t> sensor_data,
                                 int lane, double ref_vel, trajectory_t prev_path, int total_path_size,
                                 trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  /*
  Given a possible next state, generate the appropriate trajectory to realize the next state.
  */
  bool too_close = false;
  int prev_path_size = prev_path.x_vec.size();

  // For the trajectory of every sensed vehicle
  for (int i = 0; i < sensor_data.size(); i++)
  {
    // use the next point for the car sensed to deal with lag
    vehicle_t check_car;
    check_car.x = sensor_data[i].x;
    check_car.y = sensor_data[i].y;
    check_car.s = sensor_data[i].s;
    check_car.d = sensor_data[i].d;
    check_car.vx = sensor_data[i].vx;
    check_car.vy = sensor_data[i].vy;
    check_car.speed = sensor_data[i].speed;
    check_car.yaw = sensor_data[i].yaw;

    // if the sensed car is in front of us
    if (check_car.d < (2+4*lane+2) && check_car.d > (2+4*lane-2))
    {
      check_car.s += (double) prev_path_size * 0.02 * check_car.speed;

      // if the sensed car is less than 30 meters away
      if ((check_car.s > vehicle.s) && ((check_car.s - vehicle.s) < 30))
      {
        // lower velocity so we don't crash, maybe change lanes, set flags
        switch (state)
        {
          case LANE_CHANGE_LEFT:
          lane--;
          break;

          case LANE_CHANGE_RIGHT:
          lane++;
          break;

          default:
            too_close = true;
        }
      }
    }
  }

  // If the car in front was deemed too close, slow down, otherwise if we aren't near the speed limit, speed up
  if (too_close)
  {
    ref_vel -= .224; // ~5 m/s^2
  }
  else if (ref_vel < 49.5)
  {
    ref_vel += .224;
  }

  vector<double> pts_x;
  vector<double> pts_y;

  // Define our starting reference as the position of the car currently
  double ref_x = vehicle.x;
  double ref_y = vehicle.y;
  double ref_yaw = deg2rad(vehicle.yaw);

  // Start our path by looking at the previous point traversed and the current one
  if (prev_path_size < 2)
  {
    double prev_car_x = vehicle.x - cos(vehicle.yaw);
    double prev_car_y = vehicle.y - sin(vehicle.yaw);

    pts_x.push_back(prev_car_x);
    pts_x.push_back(vehicle.x);

    pts_y.push_back(prev_car_y);
    pts_y.push_back(vehicle.y);
  }
  else
  {
    ref_x = prev_path.x_vec[prev_path_size-1];
    ref_y = prev_path.y_vec[prev_path_size-1];

    double ref_x_prev = prev_path.x_vec[prev_path_size-2];
    double ref_y_prev = prev_path.y_vec[prev_path_size-2];
    ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

    pts_x.push_back(ref_x_prev);
    pts_x.push_back(ref_x);

    pts_y.push_back(ref_y_prev);
    pts_y.push_back(ref_y);
  }

  // Add even 30 meter points in front of the car's reference position to our list of next points
  vector<double> next_wp0 = getXY(vehicle.s+30, (2+lane*4), map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
  vector<double> next_wp1 = getXY(vehicle.s+60, (2+lane*4), map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
  vector<double> next_wp2 = getXY(vehicle.s+90, (2+lane*4), map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);

  pts_x.push_back(next_wp0[0]);
  pts_x.push_back(next_wp1[0]);
  pts_x.push_back(next_wp2[0]);

  pts_y.push_back(next_wp0[1]);
  pts_y.push_back(next_wp1[1]);
  pts_y.push_back(next_wp2[1]);

  // Define car's position as (0,0) with yaw 0
  for (int i = 0; i < pts_x.size(); i++)
  {
    double shift_x = pts_x[i] - ref_x;
    double shift_y = pts_y[i] - ref_y;

    pts_x[i] = (shift_x * cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
    pts_y[i] = (shift_x * sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
  }

  // Create a spline curve between the points we generated
  tk::spline s;
  s.set_points(pts_x, pts_y);

  // break up spline points so that we travel at our reference velocity
  double target_x = 30.0;
  double target_y = s(target_x);
  double target_dist = sqrt(target_x*target_x+target_y*target_y);

  double x_add_on = 0;

  trajectory_t trajectory;
  trajectory.ref_vel = ref_vel;
  trajectory.final_lane = lane;

  // Interpolate between the points we generated previously using the spline function until we have a total of 50 points
  for (int i = 1; i <= total_path_size-prev_path_size; i++)
  {
    double N = target_dist/(0.02*ref_vel/2.24);
    double x_point = x_add_on+target_x/N;
    double y_point = s(x_point);

    x_add_on = x_point;

    double x_ref = x_point;
    double y_ref = y_point;

    // Transfer car coordinates back to normal instead of relative to car
    x_point = (x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw));
    y_point = (x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw));

    x_point += ref_x;
    y_point += ref_y;

    trajectory.x_vec.push_back(x_point);
    trajectory.y_vec.push_back(y_point);
  }

  return trajectory;
}

double calculate_cost(trajectory_t trajectory, vector<trajectory_t> observed_trajectories, int lanes_available, double speed_limit,
                      trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  // Calculate cost of a trajectory by checking if:
  // - there will be a collision
  // - the car goes offroad
  // - the lane change will result in a slower velocity

  double cost = DBL_MAX;

  // if the car goes offroad or goes into the lane going in the wrong direction
  if (trajectory.final_lane < 0 || trajectory.final_lane >= lanes_available)
  {
    return cost;
  }

  // Check for trajectory collisions here
  if (trajectory.x_vec.size() > 1)
  {
    // linear estimation of cars trajectory (for speed of computation)
    point_t p1;
    p1.x = trajectory.x_vec[0];
    p1.y = trajectory.y_vec[0];

    point_t q1;
    q1.x = trajectory.x_vec[trajectory.x_vec.size()-1];
    q1.y = trajectory.y_vec[trajectory.y_vec.size()-1];

    vector<double> p1_frenet = getFrenet(p1.x, p1.y, atan2(p1.y, p1.x), map_waypoints.x_vec, map_waypoints.y_vec);
    vector<double> q1_frenet = getFrenet(q1.x, q1.y, atan2(q1.y, q1.x), map_waypoints.x_vec, map_waypoints.y_vec);

    point_t p1_left_edge;
    vector<double> p1_left_edge_temp;
    p1_left_edge_temp = getXY(p1_frenet[0], p1_frenet[1] - 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    p1_left_edge.x = p1_left_edge_temp[0];
    p1_left_edge.y = p1_left_edge_temp[1];

    point_t q1_left_edge;
    vector<double> q1_left_edge_temp;
    q1_left_edge_temp = getXY(q1_frenet[0], q1_frenet[1] - 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    q1_left_edge.x = q1_left_edge_temp[0];
    q1_left_edge.y = q1_left_edge_temp[1];

    point_t p1_right_edge;
    vector<double> p1_right_edge_temp;
    p1_right_edge_temp = getXY(p1_frenet[0], p1_frenet[1] + 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    p1_right_edge.x = p1_right_edge_temp[0];
    p1_right_edge.y = p1_right_edge_temp[1];

    point_t q1_right_edge;
    vector<double> q1_right_edge_temp;
    q1_right_edge_temp = getXY(q1_frenet[0], q1_frenet[1] + 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
    q1_right_edge.x = q1_right_edge_temp[0];
    q1_right_edge.y = q1_right_edge_temp[1];

    vector<point_t> car_edge_trajectories = { p1_left_edge, q1_left_edge, p1_right_edge, q1_right_edge };

    // account for size of cars
    for (int i = 0; i < observed_trajectories.size(); ++i)
    {

      for (int j = 0; j < observed_trajectories[i].x_vec.size()-1; ++j)
      {
        point_t p2;
        p2.x = observed_trajectories[i].x_vec[j];
        p2.y = observed_trajectories[i].y_vec[j];

        point_t q2;
        q2.x = observed_trajectories[i].x_vec[j+1];
        q2.y = observed_trajectories[i].y_vec[j+1];

        vector<double> p2_frenet = getFrenet(p2.x, p2.y, atan2(p2.y, p2.x), map_waypoints.x_vec, map_waypoints.y_vec);
        vector<double> q2_frenet = getFrenet(q2.x, q2.y, atan2(q2.y, q2.x), map_waypoints.x_vec, map_waypoints.y_vec);

        point_t p2_left_edge;
        vector<double> p2_left_edge_temp;
        p2_left_edge_temp = getXY(p2_frenet[0], p2_frenet[1] - 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
        p2_left_edge.x = p2_left_edge_temp[0];
        p2_left_edge.y = p2_left_edge_temp[1];

        point_t q2_left_edge;
        vector<double> q2_left_edge_temp;
        q2_left_edge_temp = getXY(q2_frenet[0], q2_frenet[1] - 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
        q2_left_edge.x = q2_left_edge_temp[0];
        q2_left_edge.y = q2_left_edge_temp[1];

        point_t p2_right_edge;
        vector<double> p2_right_edge_temp;
        p2_right_edge_temp = getXY(p2_frenet[0], p2_frenet[1] + 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
        p2_right_edge.x = p2_right_edge_temp[0];
        p2_right_edge.y = p2_right_edge_temp[1];

        point_t q2_right_edge;
        vector<double> q2_right_edge_temp;
        q2_right_edge_temp = getXY(q2_frenet[0], q2_frenet[1] + 1.5, map_waypoints_s, map_waypoints.x_vec, map_waypoints.y_vec);
        q2_right_edge.x = q2_right_edge_temp[0];
        q2_right_edge.y = q2_right_edge_temp[1];

        vector<point_t> obs_edge_trajectories = { p2_left_edge, q2_left_edge, p2_right_edge, q2_right_edge };

        // if there will be a collision with an observed vehicle, return the cost
        for (int k = 0; k < car_edge_trajectories.size(); k += 2) {
          for (int l = 0; l < obs_edge_trajectories.size(); l += 2) {
            if (doIntersect(car_edge_trajectories[k], car_edge_trajectories[k+1], obs_edge_trajectories[l], obs_edge_trajectories[l+1]))
            {
              // printf("traj %d, iter %d (%lf, %lf) (%lf, %lf) intersects (%lf, %lf) (%lf, %lf)\n", i, j, car_edge_trajectories[k].x,
              //                                                                                           car_edge_trajectories[k].y,
              //                                                                                           car_edge_trajectories[k+1].x,
              //                                                                                           car_edge_trajectories[k+1].y,
              //                                                                                           obs_edge_trajectories[l].x,
              //                                                                                           obs_edge_trajectories[l].y,
              //                                                                                           obs_edge_trajectories[l+1].x,
              //                                                                                           obs_edge_trajectories[l+1].y);
              return cost;
            }
          }
        }
      }
    }
  }

  // Rank trajectory based on distance from speed limit
  if (trajectory.x_vec.size() > 1)
  {
    double vx_final = trajectory.x_vec[trajectory.x_vec.size()-1] - trajectory.x_vec[trajectory.x_vec.size()-2];
    double vy_final = trajectory.y_vec[trajectory.y_vec.size()-1] - trajectory.y_vec[trajectory.y_vec.size()-2];
    double speed_final = sqrt(pow(vx_final, 2)+pow(vy_final, 2)) * 50 * 2.24;
    cost = abs(speed_limit - speed_final);
  }
  else{
    cost = speed_limit;
  }

  return cost;
}

trajectory_t choose_next_trajectory(vehicle_t vehicle, vehicle_state_t & current_state, int & lane, int lanes_available, vector<vehicle_t> sensor_data,
                                    double & ref_vel, trajectory_t prev_path, int total_path_size,
                                    trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  vector<vehicle_state_t> states = successor_states(current_state, lane, lanes_available);
  float cost;
  vector<float> costs;
  vector<trajectory_t> final_trajectories;
  vector<trajectory_t> observed_trajectories;

  for (int i = 0; i < sensor_data.size(); ++i)
  {
    observed_trajectories.push_back(constant_speed_trajectory(sensor_data[i], 50, map_waypoints, map_waypoints_del, map_waypoints_s));
  }

  for (int i = 0; i < states.size(); ++i)
  {
    trajectory_t trajectory = generate_trajectory(vehicle, states[i], sensor_data, lane, ref_vel, prev_path, total_path_size, map_waypoints, map_waypoints_del, map_waypoints_s);
    if (trajectory.x_vec.size() != 0) {
        cost = calculate_cost(trajectory, observed_trajectories, lanes_available, 50.0, map_waypoints, map_waypoints_del, map_waypoints_s);
        costs.push_back(cost);
        final_trajectories.push_back(trajectory);
    }
  }

  vector<float>::iterator best_cost = min_element(begin(costs), end(costs));
  int best_idx = distance(begin(costs), best_cost);

  current_state = states[best_idx];
  trajectory_t final_trajectory = final_trajectories[best_idx];
  ref_vel = final_trajectory.ref_vel;
  lane = final_trajectory.final_lane;
  // cout << current_state << ", " << ref_vel << ", " << lane << ", " << costs[best_idx] << endl;
  return final_trajectory;
}

/*****************************************************************************
 * Main Entrypoint
 *****************************************************************************/

int main()
{
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  trajectory_t map_waypoints;
  trajectory_t map_waypoints_del;
  vector<double> map_waypoints_s;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints.x_vec.push_back(x);
    map_waypoints.y_vec.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_del.x_vec.push_back(d_x);
    map_waypoints_del.y_vec.push_back(d_y);
  }

  // start in lane 1
  int lane = 1;

  int num_lanes = 3;

  // reference velocity to target
  double ref_vel = 0.0; // mph

  // starting state
  vehicle_state_t curr_state = KEEP_LANE;

  h.onMessage([&curr_state, &ref_vel,&map_waypoints,&map_waypoints_s,&map_waypoints_del,&lane,&num_lanes]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    // auto sdata = string(data).substr(0, length);
    // cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry")
        {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];
          vector<double> previous_path_x = j[1]["previous_path_x"];
          vector<double> previous_path_y = j[1]["previous_path_y"];
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];
          vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          /**********************************************************************
           * Build Car Trajectory
           *********************************************************************/
          trajectory_t previous_path;
          previous_path.x_vec = previous_path_x;
          previous_path.y_vec = previous_path_y;

          int prev_size = previous_path_x.size(); // number of points from the previous path which weren't traversed yet

          // If there was a previous path calculated, start calculating from the last point on that path
          if (prev_size > 0)
          {
            car_s = end_path_s;
            car_d = end_path_d;
          }

          vehicle_t this_car;
          this_car.x = car_x;
          this_car.y = car_y;
          this_car.s = car_s;
          this_car.d = car_d;
          this_car.yaw = car_yaw;
          this_car.speed = car_speed;

          // Look at sensor fusion data about the cars around us and predict their trajectories
          vector<vehicle_t> sensor_data;
          for (int i = 0; i < sensor_fusion.size(); i++)
          {
            vehicle_t check_car;
            check_car.id = sensor_fusion[i][0];
            check_car.x = sensor_fusion[i][1];
            check_car.y = sensor_fusion[i][2];
            check_car.vx = sensor_fusion[i][3];
            check_car.vy = sensor_fusion[i][4];
            check_car.s = sensor_fusion[i][5];
            check_car.d = sensor_fusion[i][6];
            check_car.speed = sqrt(pow(check_car.vx, 2)+pow(check_car.vy, 2));
            check_car.yaw = atan2(check_car.vy, check_car.vx);

            // Generate a trajectory assuming contant speed for all cars detected through sensor fusion
            sensor_data.push_back(check_car);
          }

          trajectory_t next_trajectory;

          // Add untraversed points calculated in the previous iteration to this path
          for (int i = 0; i < previous_path_x.size(); i++)
          {
            next_trajectory.x_vec.push_back(previous_path_x[i]);
            next_trajectory.y_vec.push_back(previous_path_y[i]);
          }

          // Choose next trajectory for our car
          trajectory_t new_trajectory;
          new_trajectory = choose_next_trajectory(this_car, curr_state, lane, num_lanes, sensor_data, ref_vel, previous_path, 50, map_waypoints, map_waypoints_del, map_waypoints_s);

          // Add newly calculated trajectory to this path
          merge_trajectories(next_trajectory, new_trajectory);

          // Send New Trajectory
          json msgJson;
          msgJson["next_x"] = next_trajectory.x_vec;
          msgJson["next_y"] = next_trajectory.y_vec;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      }
      else
      {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
