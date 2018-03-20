# Implenting the Path Planner for the Udacity Term 3 Simulator

### By Steven Eisinger

## General Overview

The goal of this project is to plan the path of a simulated vehicle along a highway without crashing, breaking any laws, or causing passenger discomfort. The main program flow goes as follows:

1. Vehicle observations, map data, and the currently followed trajectory are passed in and coverted into objects.
1. The trajectory passed in to the simulator is added to our newly computed trajectory.
1. In `choose_next_trajectory`: Next possible states are checked.
1. In `choose_next_trajectory`: Each observed vehicle has a constant speed trajectory generated with `constant_speed_trajectory`.
1. In `choose_next_trajectory`: A new trajectory for each possible next state is computed for the number of points missing from the passed in trajectory with `generate_trajectory`.
1. In `choose_next_trajectory`: Each trajectory is given a cost with `calculate_cost`.
1. The trajectory with the lowest cost is returned and appended to the original trajectory that was passed in.
1. The new trajectory is sent to the simulator and the program waits for further input.
1. When data are received, begin at step 1.

The sections below will go over the functions described above in detail. I will go in an order which improves understanding.

## Generating a Trajectory

Trajectory calculation could be done in numerous different ways, such as generating many possible trajectories for a probabilistic end location and choosing the best one or staying in one lane and checking that the car never crashes or goes over the speed limit. This program uses the spline technique as suggested by the project walkthrough. Nearly all of this happens in the `generate_trajectory` function. It goes as follows:

1. We start the new trajectory from the last point of the previous trajectory passed into the program.
1. We pass in a next possible state for the vehicle.
1. Every car that was detected is checked. If there is a car in from of us, we change the intended final lane for the trajectory or slow down.
    ```c++
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
    ```

1. Adjust the car's current speed based on state.
    ```c++
        // If the car in front was deemed too close, slow down, otherwise if we aren't near the speed limit, speed up
        if (too_close)
        {
            ref_vel -= .224; // ~5 m/s^2
        }
        else if (ref_vel < 49.5)
        {
            ref_vel += .224;
        }
    ```

1. Generate points for a spline based on the car's position at the end of the previouslt computed trajectory and 3 30 meter intervals further is s (in Frenet coordinates).
    ```c++
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
    ```

1. Transform the points so that the car is at (0, 0).
    ```c++
        // Define car's position as (0,0) with yaw 0
        for (int i = 0; i < pts_x.size(); i++)
        {
            double shift_x = pts_x[i] - ref_x;
            double shift_y = pts_y[i] - ref_y;

            pts_x[i] = (shift_x * cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
            pts_y[i] = (shift_x * sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
        }
    ```

1. Create a spline and interpolate through the spline such that the car will travel at the desired speed.
    ```c++
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
    ```

1. Return the trajectory for further computation.

This is the heart of the program. The spline guarantees smooth actuations through the trajectory and makes computation simple because the spline function has been imported from _spline.h_ for us.

There is a similar function, called `constant_speed_trajectory`, which does something extremely similar for all cars observed by our vehicle, with the following differences:

1. There is only one state: KEEP_LANE.
1. The car might be going in the opposite direction, and has to be accounted for.
1. A car might be seen as not moving, and its trajectory should accomodate that.
1. The car is assumed to go at a constant speed.
1. A full 50 point trajectory is always computed, so that we don't need to keep track of the trajectories of all of these cars (but this requires more computation time from the program).

`constant_speed_trajectory` is called on all observed vehicles so that it can be used in combination with our vehicle's trajectory to compute a cost.

## Computing Trajectory Cost

