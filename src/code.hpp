//
//  code.cpp
//  path_planning
//
//  Created by LI JIA on 20.04.19.
//

#ifndef PATHPLANNING_CODE_H_
#define PATHPLANNING_CODE_H_
#include "helpers.h"
#include "spline.h"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include <cassert>


namespace pp {
  using namespace std;
  using namespace Eigen;

  static const int frame_rate = 1000/20; //< 50 fps
  static const double frame_time = 0.020; //< 20ms
  static const double speed_limit = 22.35; //< m/s
  static const double speed_limit_per_frame = speed_limit / frame_rate;
  static const double accel_limit = 10.; //< m/s^2
  static const int max_output_size = frame_rate * 1; //< output 1 second of waypoints
  static const double max_yaw_rate = 0.;

  double target_lane_ = 0.;
  int cur_frame = 0;
  int last_commit_size = 0;

  class PathGenerator {
    const vector<double>& map_x_;
    const vector<double>& map_y_;
    const vector<double>& map_s_;
    const vector<double>& map_dx_;
    const vector<double>& map_dy_;

  public:
    PathGenerator(
        const vector<double>&  map_waypoints_x,
        const vector<double>&  map_waypoints_y,
        const vector<double>&  map_waypoints_s,
        const vector<double>&  map_waypoints_dx,
        const vector<double>&  map_waypoints_dy)
      : map_x_(map_waypoints_x),
      map_y_(map_waypoints_y),
      map_s_(map_waypoints_s),
      map_dx_(map_waypoints_dx),
      map_dy_(map_waypoints_dy)
    {
    }

    void GenPath(Vector2d car_pt, double s_car, double d_car, Vector2d prev_target_fernet, double yaw_deg, double speed, vector<double>& out_xs,
        vector<double>& out_ys) {
      double yaw = deg2rad(yaw_deg);
      
      auto n_committed = out_xs.size();
      assert(out_ys.size() == n_committed);
      
      // generate at least 3 frames a time
      if (max_output_size - n_committed < 3)
        return;
      
      cur_frame += last_commit_size - n_committed;

      // continue from last committed point
      if (n_committed > 0)
      {
        size_t idx_last = n_committed;
        car_pt = Vector2d(out_xs[idx_last], out_ys[idx_last]);
        yaw = atan2(out_ys[idx_last]-out_ys[idx_last-1], out_xs[idx_last]-out_xs[idx_last-1]);
      }
      else
      {
        prev_target_fernet = Vector2d{s_car, d_car};
      }

      assert(speed >= 0.);
      //int wpidx_a = ClosestWaypoint(car_pt.x(), car_pt.y(), map_x_, map_y_);
      //int wpidx_b = GetPrevWp(car_fernet.x());
      //assert( abs(wpidx_a - wpidx_b) <= 1);

      double yaw_rate = pi() / 2.;   //< 1/4 * 2pi
      double target_accel = accel_limit;
      double target_yaw = yaw + frame_time * yaw_rate;
      double next_speed = min(speed_limit, speed + frame_time * target_accel);
      double d = 4. * target_lane_ + 2.;

      // === step x: generate jerk minimizing trajectory ===
      // input: 

      auto T = frame_time * (max_output_size - n_committed);
      auto dt_commit = frame_time * n_committed;
      auto s_init_p = prev_target_fernet.x();
      auto s_init_v = dt_commit > 0 ? (s_init_p - s_car) / dt_commit : 0.0;
      auto s_init_a = dt_commit > 0 ? (s_init_v - .9*speed) / dt_commit : 0.0;

      // accelerate to start
      auto s_final_a = (s_init_v <= .5 * speed_limit) ? accel_limit : 0.0;

      auto s_final_v = min(speed_limit, s_init_v + s_final_a * T);
      auto s_final_p = s_init_p + s_init_v*T +.5*s_final_a*T*T;

      auto args_s = JMT(vector<double>{s_init_p, s_init_v, s_init_a}, vector<double>{s_final_p, s_final_v, s_final_a}, T);

      vector<double> d_init_vec={prev_target_fernet.y(), 0, 0}, d_final_vec={4.0 * target_lane_ + 2.0, 0, 0};
      auto args_d = JMT(d_init_vec, d_final_vec, T);

      // get spline-based fernet to xy projection
      //auto sp_pair = GetSpline(s_init_p, s_final_p, d_init_vec[0], d_final_vec[0]);

      for (int i = 0; i < (max_output_size - n_committed); ++i)
      {
        auto s = EvalPoly(args_s, (i + 1) * frame_time);
        auto d = EvalPoly(args_d, (i + 1) * frame_time);

        auto xy = ::getXY(s,d,map_s_, map_x_,map_y_);
        auto x = xy[0], y = xy[1];
        // double x = sp_pair.first(s), y = sp_pair.second(s);
        out_xs.push_back(x);
        out_ys.push_back(y);
      }

      if (cur_frame > 300)
        CheckOutput(out_xs, out_ys);

      last_commit_size = out_xs.size();

      // auto sp = GetSpline(s_start, d);

      // for (int i=1; i < max_output_size - n_committed; ++i)
      // {
      //   // next_pt = car_pt + target_speed * dt * Vector2d(cos(yaw+ dt* yaw_rate), sin(yaw+dt*yaw_rate));
      //   Vector2d next_fernet = target_fernet + Vector2d(0.5 * (i+1),0);//next_speed * dt * (i+1), 0);
      //   double next_s  = next_fernet.x();

      //   // auto xy = ::getXY(next_fernet.x(),next_fernet.y(),map_s_, map_x_,map_y_);
      //   double next_x = sp.first(next_s);
      //   double next_y = sp.second(next_s);
      //   out_xs.push_back(next_x);
      //   out_ys.push_back(next_y);

      //   // update variables
      //   // speed = target_speed;
      //   // yaw = target_yaw;
      //   // // car_pt = next_pt;
      //   // car_fernet = next_fernet;
      // }
    }

    void CheckOutput(const vector<double>& out_xs, const vector<double>& out_ys)
    {
      double dt = 0.02;
      for (int i =0; i< out_xs.size()-2;++i)
      {
        auto x0 = out_xs[i];
        auto y0 = out_ys[i];
        auto x1 = out_xs[i+1];
        auto y1 = out_ys[i+1];
        auto x2 = out_xs[i+2];
        auto y2 = out_ys[i+2];

        auto vel0 = distance(x0,y0,x1,y1) / dt;
        auto vel1 = distance(x1,y1,x2,y2) / dt;

        auto acc = (vel1 - vel0) / dt;

        assert(vel0 < 100 && vel1 < 100);
        assert(abs(acc) < 100);
      }
    }


    int GetPrevWp(double s)
    {
      int prev_wp = -1;

      while (s > map_s_[prev_wp+1] && (prev_wp < (int)(map_s_.size()-1))) {
        ++prev_wp;
      }
      return prev_wp;
    }

    // Vector2d GetCoords(Vector2d fernet)
    // {
    //   auto sp = GetSpline(fernet.x(), fernet.y());
    //   return {sp.first(fernet.x()), sp.second(fernet.x())};
    // }

    pair<tk::spline, tk::spline> GetSpline(double s_start, double s_end, double d_start, double d_end)
    {
      int idx0 = GetPrevWp(s_start);

      int idx_n1 = (idx0 - 1 + map_s_.size()) % map_s_.size();
      int idx_n2 = (idx0 - 2 + map_s_.size()) % map_s_.size();
      int idx_n3 = (idx0 - 3 + map_s_.size()) % map_s_.size();

      int idx1 = (idx0 + 1) % map_s_.size();
      int idx2 = (idx0 + 2) % map_s_.size();
      int idx3 = (idx0 + 3) % map_s_.size();

      int indices[7] = {idx_n3, idx_n2, idx_n1, idx0, idx1, idx2, idx3};

      vector<double> xs, ys, inputs;
      for (int i = 1; i < 6; ++i)
      {
        int idx = indices[i];
        double s = map_s_[idx];
        double d;
        if (s <= s_start)
          d = d_start;
        else if (s < s_end)
          d = d_start + (((s - s_start)/(s_end - s_start)) * (d_end - d_start));
        else
          d = d_end;

        auto estimate_xy = ::getXY(s, d, map_s_, map_x_, map_y_);
        xs.push_back(estimate_xy[0]);
        ys.push_back(estimate_xy[1]);

        double x_prev = map_x_[indices[i - 1]];
        double y_prev = map_y_[indices[i - 1]];
        double x_next = map_x_[indices[i + 1]];
        double y_next = map_y_[indices[i + 1]];
        double x_wp = map_x_[idx], y_wp = map_y_[idx];

        double heading_prev = atan2(y_wp - y_prev, x_wp - x_prev);
        double heading = atan2(y_next - y_wp, x_next - x_wp);
        double heading_avg = (heading_prev + heading) / 2.;
        double prep = heading_avg - pi() / 2.;
        double x = x_wp + d * cos(prep);
        double y = y_wp + d * sin(prep);

        // assert(estimate_xy[0] - x < 1. && estimate_xy[1] - y < 1.);

        // xs.push_back(x);
        // ys.push_back(y);
        inputs.push_back(s);
      }

      tk::spline sp_x, sp_y;
      sp_x.set_points(inputs, xs);
      sp_y.set_points(inputs, ys);

      return {sp_x,sp_y};
    }


  private:
    static vector<double> JMT(const vector<double> &start, const vector<double> &end, double T)
    {
      Eigen::MatrixXd m(3,3);

      m <<
        pow(T,3), pow(T,4), pow(T,5),
        3.*pow(T,2), 4.*pow(T,3), 5.*pow(T,4),
        6.*T, 12.*pow(T,2), 20.*pow(T,3) ;

      m = m.inverse();

      auto rhs = Vector3d{
          end[0] - (start[0] + start[1] * T + .5 * start[2] * T * T),
          end[1] - (start[1] + start[2] * T),
          end[2] - start[2]};

      auto br = m * rhs;

      return {start[0],start[1],.5*start[2], br[0],br[1],br[2]};
    }

    static double EvalPoly(const vector<double>& a, double t)
    {
      return a[0] + a[1] * t + a[2] * t * t + a[3] * t * t * t + a[4] * pow(t, 4) + a[5] * pow(t, 5);
    }
  };


}

#endif // PATHPLANNING_CODE_H_
// vim: sw=2 ts=2 sts=2 et:
