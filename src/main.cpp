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

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

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
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
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

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
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
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
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

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

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
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // 0,1,2 == left,middle,right lanes; start in lane 1
  int lane = 1;
  // Move a reference velocity to target
  double ref_vel = 0.0; // mph

  // h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
  //                    uWS::OpCode opCode) {
  // h.onMessage([&ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&lane_change_wp](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
  //                    uWS::OpCode opCode) {
  h.onMessage([&ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          // ego's localization data
        	double car_x = j[1]["x"];
        	double car_y = j[1]["y"];
        	double car_s = j[1]["s"];
        	double car_d = j[1]["d"];
        	double car_yaw = j[1]["yaw"];
        	double car_speed = j[1]["speed"];

        	// Previous path data given to the Planner
        	auto previous_path_x = j[1]["previous_path_x"];
        	auto previous_path_y = j[1]["previous_path_y"];
        	// Previous path's end s and d values 
        	double end_path_s = j[1]["end_path_s"];
        	double end_path_d = j[1]["end_path_d"];

        	// Sensor Fusion Data, a list of all other cars on the same side of the road.
        	auto sensor_fusion = j[1]["sensor_fusion"];
          //vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          int prev_size = previous_path_x.size();

          //beginning of commented-out code from demo
          // sensor fusion
          if (prev_size > 0) {
            car_s = end_path_s;
          }

          bool car_directly_ahead_too_close = false;
          bool car_tailgating = false;
          bool option_to_move_left = true;
          bool option_to_move_right = true;

          //std::cout << "Number of other cars on the road:" << sensor_fusion.size();
          // find ref_v to use
          // for all cars in the sensor fusion data
          for (int i = 0; i < sensor_fusion.size(); i++) {
            //std::cout << "i = " << i;
            // i corresponds to ith car on the road
            float d = sensor_fusion[i][6]; // horizontal position of ith car
            double vx = sensor_fusion[i][3];
            double vy = sensor_fusion[i][4];
            double check_speed = sqrt(vx*vx+vy*vy);
            double check_car_s = sensor_fusion[i][5]; // vertical position of ith car

            // magic numbers
            // There are formulas for calculating a safe distance between vehicles at highway speed.
            // They involve the speed and lengths of vehicles
            // Assuming a speed of 50 mph, which is, 73 ft/sec, and a vehicle length of 16 feet,
            // the formula would have a driver leave 1 second for every ten feet in vehicle length plus one second,
            // resulting in: (1.6 + 1) seconds * 73.3 ft/s = 191 ft = 58 m
            double min_safe_dist = 20.0;
            // double trailing_min_safe_dist = 20.0;

            // L3M8
            // With Frenet coordinates, we use the variables x and d to describe a vehicle's position on the road.
            // The s coordinate represents distance along the road (also known as longitudinal displacement).
            // The d coordinate represents side-to-side position on the road (also known as lateral displacement).
            // A typical trajectory in Frenet coordinates looks like a straight line.
            // If the vehicle were moving at a constant speed, we could write a mathematical description of the 
            // vehicle's position.
            // Straight lines are much easier to work with; use Frenet.
            // d represents horizontal distance or lane
            // s represents vertical distance or ahead or behind ego

            check_car_s += (double)prev_size*0.02*check_speed; // if using previous points can projet s value out

            // Loop through all of the other cars to determine if they are in ego's lane or too close to ego or blocking ego
            // from changing lanes.
            //
            // lane width is 4m, therefore, 1/2 lane width is 2m to get lane center
            // lane is 0, 1, or 2 for left, middle, right lanes
            // 0 < d < 4; lane 0
            // 4 < d < 8; lane 1
            // 8 < d < 12; lane 2

            // lane edges
            double lane_left_edge = 2+4*lane-2;
            double lane_right_edge = 2+4*lane+2;

            // if the ith car is in the same lane as ego...
            if (d > lane_left_edge && d < lane_right_edge) {
              // if ith car is ahead of ego...
              if (check_car_s > car_s) {
                // if the distance is not safe...
                if ((check_car_s - car_s) < min_safe_dist) {
                  // will eventually slow down or change lanes...
                  // TODO: fix this using finite state machine and frenet
                  car_directly_ahead_too_close = true;
                } // car too near ego
              } // car ahead of ego
              else if (check_car_s < car_s) {
                // if the distance is not safe...
                if ((car_s - check_car_s) < min_safe_dist) {
                  // will eventually slow down or change lanes...
                  // TODO: fix this using finite state machine and frenet
                  car_tailgating = true;
                } // car too near ego
              } // car ahead of ego
            } // car in same lane as ego

            // Check if either the lane to the left or right of ego is clear to move into.

            // investigate changing lanes
            double s_distance = abs(car_s-check_car_s);
            bool car_faster_than_ego = check_speed > car_speed;

            // investigate moving left
            
            // if there is a lane to the left of ego...
            if (lane > 0) {
              // lane to the left of ego
              int lane_left_of_ego = lane-1;
              // to-left lane edges
              double to_left_lane_left_edge = 2+4*lane_left_of_ego-2;
              double to_left_lane_right_edge = 2+4*lane_left_of_ego+2;

              // process cars to the left
              // if the ith car is in the lane to the left of ego...
              if (d > to_left_lane_left_edge && d < to_left_lane_right_edge)
                // if ((check_speed > car_speed) && (check_car_s < car_s)) {
                //  min_safe_distance_behind = 50;
                //  std::cout << "car approaching on left" << std::endl;
                // }

                // The lane is not available if the car in the lane to the left of ego is going faster than ego or if it is
                // within a non-safe s distance of ego.

                // if car to the left is not within a safe distance or it is going faster than ego...
                if (s_distance < min_safe_dist || car_faster_than_ego) {
                  option_to_move_left = false;
                } // is lane to the left available
            } // a lane to the left exists

            // investigate moving right

            // if there is a lane to the right of ego...
            if (lane < 2) {
              // lane to the right of ego
              int lane_right_of_ego = lane+1;
              // to-right lane edges
              double to_right_lane_left_edge = 2+4*lane_right_of_ego-2;
              double to_right_lane_right_edge = 2+4*lane_right_of_ego+2;

              // process cars to the right
              // if the ith car is in the lane to the right of ego...
              if (d > to_right_lane_left_edge && d < to_right_lane_right_edge)
                // if ((check_speed > car_speed) && (check_car_s < car_s)) {
                //  min_safe_distance_behind = 50;
                //  std::cout << "car approaching on left" << std::endl;
                // }

                // The lane is not available if the car in the lane to the right of ego is going faster than ego or if it is
                // within a non-safe s distance of ego.

                // if car to the right is not within a safe distance or it is going faster than ego...
                if (s_distance < min_safe_dist || car_faster_than_ego) {
                  option_to_move_right = false;
                } // is lane to the right available
            } // a lane to the right exists
          } // loop through other cars

          /* At this point, we know the positions of the cars ane if we should:
           * remain in our lane at same speed,
           * remain in our lane and alter speed (slower or faster),
           * try to change lanes
           */

          // if ith car is directly ahead of ego and too close...
          if (car_directly_ahead_too_close) {
            std::cout << "Car directly ahead and too close!" << std::endl;
            // 1 mi = 1609.34 m
            // ref_vel -= 0.224; // creates a deceleration of approx 5/ms2 (converted from m/s to mph)

            // See if a lane is available for ego to change to.

            // print options to console
            // if there exists a lane to the left of ego and it is available...
            if (lane > 0 && option_to_move_left) {
              std::cout << "  option to move left" << std::endl;
            } else {
              std::cout << "  not an option to move left" << std::endl;
            }

            // if there exists a lane to the right of ego and it is available
            if (lane < 2 && option_to_move_right) {
              std::cout << "  option to move right" << std::endl;
            } else {
              std::cout << "  not an option to move right" << std::endl;
            }

            // if there exists a lane to the left of ego and it is available...
            if (lane > 0 && option_to_move_left) {
              std::cout << "    ego moving left..." << std::endl;
              lane -= 1;
            } else if (lane < 2 && option_to_move_right) {
              std::cout << "    ego moving right..." << std::endl;
              lane += 1;
            } else {
              std::cout << "    ego is stuck in this lane and is decelerating" << std::endl;
              ref_vel -= 0.1;
            }
          }
          else if (ref_vel < 49.5) {
            std::cout << "open road ahead; ego is accelerating" << std::endl;
            ref_vel += 0.224;
          }
          //end of commented-out code

        	json msgJson;

        	vector<double> next_x_vals;
        	vector<double> next_y_vals;


          // This starts from the classroom lesson. The initial version of the code moves the car forward.
          // double dist_inc = 0.5;
          // for(int i = 0; i < 50; i++)
          // {
          //   next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad(car_yaw)));
          //   next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw)));
          // }
          // The second version of this code uses Frenet coordinates to keep the car in its lane.
        	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          // double dist_inc = 0.5;
          // for(int i = 0; i < 50; i++)
          // {
          //   // we want to do i+1; otherwise, our car will be exactly where the first point is at
          //   // and it won't be transitioning - it will be sitting still.
          //   double next_s = car_s+(i+1)*dist_inc;
          //   // we are in the middle lane
          //   // the waypoints are measured from the double yellow line in the middle of the road
          //   // so we are like 1 1/2 lanes from where the waypoints are
          //   // from the classroom, lanes are four meters wide
          //   double next_d = 6; // (1.5 lanes * 4 m/lane)
          //   vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

          //   next_x_vals.push_back(xy[0]);
          //   next_y_vals.push_back(xy[1]);    
          // }
          // double dist_inc = 0.3; // the distance between points; 0.5 (m) is close to 50 mph
          // 50 points for the path planner
          // for(int i = 0; i <50; i++)
          // {
          //   double next_s = car_s+(i+1)*dist_inc;
          //   double next_d = 6;
          //   vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

          //   next_x_vals.push_back(xy[0]);
          //   next_y_vals.push_back(xy[1]);
          // }

          // END
         //  msgJson["next_x"] = next_x_vals;
        	// msgJson["next_y"] = next_y_vals;


          // Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
          // Later we will interpolate these waypoints with a spline and fill it in with more points that control speed
          vector<double> ptsx;
          vector<double> ptsy;

          // referenece x,y, yaw states
          // either we will reference the starting point or where the car is at previous paths end point
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          // if previous size is almost empty, use the car as starting reference
          if (prev_size < 2) {
            // Use two points that make the path tangent to the car
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);

            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          }
          // use the previous path's end point as starting reference
          else {
            // Redefine reference state as previous path end point
            ref_x = previous_path_x[prev_size-1];
            ref_y = previous_path_y[prev_size-1];

            double ref_x_prev = previous_path_x[prev_size-2];
            double ref_y_prev = previous_path_y[prev_size-2];
            ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);

            // Use two points that make the path tangent to the previous path's point
            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }

          // beginning of edit
          // In Frenet add evenly 30m spaced points ahead of the starting reference
          vector<double> next_mp0 = getXY(car_s+30,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
          vector<double> next_mp1 = getXY(car_s+60,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);
          vector<double> next_mp2 = getXY(car_s+90,(2+4*lane),map_waypoints_s,map_waypoints_x,map_waypoints_y);

          ptsx.push_back(next_mp0[0]);
          ptsx.push_back(next_mp1[0]);
          ptsx.push_back(next_mp2[0]);

          ptsy.push_back(next_mp0[1]);
          ptsy.push_back(next_mp1[1]);
          ptsy.push_back(next_mp2[1]);

          // shift and rotation of the coordinate system
          for (int i = 0; i < ptsx.size(); i++ ) {
            // shift car reference angle to 0 degrees
            double shift_x = ptsx[i]-ref_x;
            double shift_y = ptsy[i]-ref_y;

            ptsx[i] = (shift_x *cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
            ptsy[i] = (shift_x *sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
          }

          // create a spline
          tk::spline s;

          // set (x,y) points to the spline
          s.set_points(ptsx,ptsy); // five anchor points

          // Define the actual (x,y) points we will use for the planner
          // vector<double> next_x_vals;
          // vector<double> next_y_vals;

          // Start with all of the previous path points from last time
          for (int i = 0; i < previous_path_x.size(); i++) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // Calculate how to break up spline points so that we travel
          double target_x = 30.0; // our horizon
          double target_y = s(target_x);
          double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

          double x_add_on = 0;

          // Fill up the rest of our path planner after filling it with previous points, here we will always output 50 points
          for (int i = 1; i <= 50-previous_path_x.size(); i++) {
            double N = (target_dist/(.02*ref_vel/2.24));
            double x_point = x_add_on+(target_x)/N;
            double y_point = s(x_point);

            x_add_on = x_point;

            double x_ref = x_point;
            double y_ref = y_point;

            // rotate back to normal after rotating it earlier
            x_point = (x_ref *cos(ref_yaw)-y_ref*sin(ref_yaw));
            y_point = (x_ref *sin(ref_yaw)+y_ref*cos(ref_yaw));
            x_point += ref_x;
            y_point += ref_y;

            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);
          }
          // end of edit

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

        	auto msg = "42[\"control\","+ msgJson.dump()+"]";

        	//this_thread::sleep_for(chrono::milliseconds(1000));
        	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
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