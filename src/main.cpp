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
    states.push_back(PREP_LANE_CHANGE_LEFT);
    states.push_back(PREP_LANE_CHANGE_RIGHT);
  }
  else if (state == PREP_LANE_CHANGE_LEFT)
  {
    if (lane != lanes_available - 1)
    {
      states.push_back(PREP_LANE_CHANGE_LEFT);
      states.push_back(LANE_CHANGE_LEFT);
    }
  }
  else if (state == PREP_LANE_CHANGE_RIGHT)
  {
    if (lane != 0)
    {
      states.push_back(PREP_LANE_CHANGE_RIGHT);
      states.push_back(LANE_CHANGE_RIGHT);
    }
  }
  //If state is LANE_CHANGE_LEFT or LANE_CHANGE_RIGHT, then just return KEEP_LANE
  return states;
}

// Used for calculating trajectories of observed vehicles
trajectory_t constant_speed_trajectory(vehicle_t vehicle)
{
  // TODO:
  trajectory_t trajectory;
  trajectory.x_vec.push_back(vehicle.x);
  trajectory.y_vec.push_back(vehicle.y);

  return trajectory;
}

trajectory_t keep_lane_trajectory(vehicle_t vehicle, vector<vehicle_t> sensor_data,
                                  int & lane, double & ref_vel, trajectory_t prev_path, int total_path_size,
                                  trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  // TODO: make this work
  bool too_close = false;
  int prev_path_size = prev_path.x_vec.size();

  // For the trajectory of every sensed vehicle
  for (int i = 0; i < sensor_data.size(); i++)
  {
    // use the next point for the car sensed to deal with lag
    vehicle_t check_car;
    vector<double> check_car_frenet = getFrenet(check_car.x, check_car.y, 0, map_waypoints.x_vec, map_waypoints.y_vec);
    check_car.x = sensor_data[i].x;
    check_car.y = sensor_data[i].y;
    check_car.s = sensor_data[i].s;
    check_car.d = sensor_data[i].d;
    check_car.vx = sensor_data[i].vx;
    check_car.vy = sensor_data[i].vy;

    // if the sensed car is in front of us
    if (check_car.d < (2+4*lane+2) && check_car.d > (2+4*lane-2))
    {
      double check_speed = sqrt(pow(check_car.vx, 2) + pow(check_car.vy, 2));
      check_car.s += (double) prev_path_size * 0.02 * check_speed;

      // if the sensed car is less than 30 meters away
      if ((check_car.s > vehicle.s) && ((check_car.s - vehicle.s) < 30))
      {
        // lower velocity so we don't crash, maybe change lanes, set flags
        too_close = true;
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

trajectory_t lane_change_trajectory(vehicle_t vehicle, vehicle_state_t state, vector<vehicle_t> sensor_data,
                                    int & lane, double & ref_vel, trajectory_t prev_path, int total_path_size,
                                    trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  // TODO:
  trajectory_t trajectory;
  return trajectory;
}

trajectory_t prep_lane_change_trajectory(vehicle_t vehicle, vehicle_state_t state, vector<vehicle_t> sensor_data,
                                         int & lane, double & ref_vel, trajectory_t prev_path, int total_path_size,
                                         trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  // TODO:
  trajectory_t trajectory;
  return trajectory;
}

trajectory_t generate_trajectory(vehicle_t vehicle, vehicle_state_t state, vector<vehicle_t> sensor_data,
                                 int & lane, double & ref_vel, trajectory_t prev_path, int total_path_size,
                                 trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  /*
  Given a possible next state, generate the appropriate trajectory to realize the next state.
  */
  trajectory_t trajectory;
  if (state == KEEP_LANE) {
      trajectory = keep_lane_trajectory(vehicle, sensor_data, lane, ref_vel, prev_path, total_path_size, map_waypoints, map_waypoints_del, map_waypoints_s);
  } else if (state == LANE_CHANGE_LEFT || state == LANE_CHANGE_RIGHT) {
      trajectory = lane_change_trajectory(vehicle, state, sensor_data, lane, ref_vel, prev_path, total_path_size, map_waypoints, map_waypoints_del, map_waypoints_s);
  } else if (state == PREP_LANE_CHANGE_LEFT || state == PREP_LANE_CHANGE_RIGHT) {
      trajectory = prep_lane_change_trajectory(vehicle, state, sensor_data, lane, ref_vel, prev_path, total_path_size, map_waypoints, map_waypoints_del, map_waypoints_s);
  }
  return trajectory;
}

double calculate_cost(trajectory_t trajectory, vector<vehicle_t> sensor_data)
{
  // TODO: calculate cost of a trajectory
  return 1.0;
}

trajectory_t choose_next_trajectory(vehicle_t vehicle, vehicle_state_t & current_state, int & lane, int lanes_available, vector<vehicle_t> sensor_data,
                                    double & ref_vel, trajectory_t prev_path, int total_path_size,
                                    trajectory_t map_waypoints, trajectory_t map_waypoints_del, vector<double> map_waypoints_s)
{
  vector<vehicle_state_t> states = successor_states(current_state, lane, lanes_available);
  float cost;
  vector<float> costs;
  vector<trajectory_t> final_trajectories;

  for (int i = 0; i < states.size(); ++i)
  {
    trajectory_t trajectory = generate_trajectory(vehicle, states[i], sensor_data, lane, ref_vel, prev_path, total_path_size, map_waypoints, map_waypoints_del, map_waypoints_s);
    if (trajectory.x_vec.size() != 0) {
        cost = calculate_cost(trajectory, sensor_data);
        costs.push_back(cost);
        final_trajectories.push_back(trajectory);
    }
  }

  vector<float>::iterator best_cost = min_element(begin(costs), end(costs));
  int best_idx = distance(begin(costs), best_cost);

  // current_state = states[best_idx];
  current_state = KEEP_LANE; // TODO: remove

  return final_trajectories[best_idx];
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
