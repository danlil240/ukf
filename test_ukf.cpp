//
// Created by ari on 8/2/22.
//

#include "UKF.hpp"
#include "model.hpp"
#include <iostream>
#include <fstream>
#include "ros_common/matplotlib-cpp/matplotlibcpp.h"
#include <random>
#include <yaml-cpp/yaml.h>
#include <ros_common/helpers/helper_methods.hpp>

namespace plt = matplotlibcpp;
bool received_gimbal_measurement_{};
bool received_camera_measurement_{};
bool received_lrf_measurement_{};

Vector8d y_{}; // Actual measurement
Vector11d x_{}; // State



std::vector<double> time_for_plot;
std::vector<double> time_for_plot_gimbal;
std::vector<double> time_for_plot_vision;
std::vector<double> time_for_plot_range;
std::vector<double> time_for_plot_gps_los;

std::vector<double> range_covariance;
std::vector<double> range_rate_covariance;
std::vector<double> pitch_covariance;
std::vector<double> pitch_rate_covariance;
std::vector<double> yaw_covariance;
std::vector<double> yaw_rate_covariance;
std::vector<double> estimated_range;
std::vector<double> estimated_range_rate;
std::vector<double> estimated_pitch;
std::vector<double> estimated_yaw;
std::vector<double> estimated_pitch_rate;
std::vector<double> estimated_yaw_rate;
std::vector<double> range_measurements;
std::vector<double> pitch_measurements;
std::vector<double> yaw_measurements;
std::vector<double> gyro_y_measurements;
std::vector<double> gyro_z_measurements;
std::vector<double> center_x_measurements;
std::vector<double> center_y_measurements;
std::vector<double> gps_los_range;
std::vector<double> gps_los_range_rate;
std::vector<double> gps_los_pitch;
std::vector<double> gps_los_pitch_rate;
std::vector<double> gps_los_yaw;
std::vector<double> gps_los_yaw_rate;
Vector11d state;
MatrixXd R(8, 8);
double pitch_measurement, yaw_measurement, center_x, center_y, range_measurement;
bool range_rate_valid = false;
bool range_valid = false;

Vector11d covariance;

double previous_lrf_measurement = 0;
double previous_lrf_time = 0;

double max_range_rate_validity;
double min_range_rate_validity;
double max_range_validity;
double min_range_validity;
std::shared_ptr<model> model_ = std::make_shared<model>();
UKF ukf;

void updateVisionMeasurement(double x1c, double y1c)
{
    
    y_(MEASUREMENT_X1C) = x1c;
    y_(MEASUREMENT_Y1C) = y1c;
    
    // Initialize LOS angles on first reading if already received gimbal measurement
    if (not received_camera_measurement_ and received_gimbal_measurement_)
    {
        model::LOSFromGimbalAndCamera(x_(STATE_PITCH), x_(STATE_YAW), y_);
    }
    
    received_camera_measurement_ = true;
    
    ukf.setMeasurementModel(model::stateToMeasurementVision);
    ukf.setMeasurementNoise(R.block(MEASUREMENT_X1C, MEASUREMENT_X1C, 2, 2));
    
    Vector11d u;
    u.setZero();
    ukf.predict(u);
    ukf.update(y_.segment(MEASUREMENT_X1C, 2));
}

void updateGimbalMeasurement(double pitch, double yaw, Vector3d gyro)
{
    y_(MEASUREMENT_GIMBAL_PITCH) = pitch;
    y_(MEASUREMENT_GIMBAL_YAW) = helpers::unwrap(y_(MEASUREMENT_GIMBAL_YAW), yaw);
    y_(MEASUREMENT_GYRO_X) = gyro.x();
    y_(MEASUREMENT_GYRO_Y) = gyro.y();
    y_(MEASUREMENT_GYRO_Z) = gyro.z();
    
    // Initialize gimbal state on first reading. Initialize LOS angles if already received camera measurement
    if (not received_gimbal_measurement_)
    {
        x_(STATE_GIMBAL_PITCH) = y_(MEASUREMENT_GIMBAL_PITCH);
        x_(STATE_GIMBAL_YAW) = y_(MEASUREMENT_GIMBAL_YAW);
        x_(STATE_GYRO_X) = y_(MEASUREMENT_GYRO_X);
        x_(STATE_GYRO_Y) = y_(MEASUREMENT_GYRO_Y);
        x_(STATE_GYRO_Z) = y_(MEASUREMENT_GYRO_Z);
        
        if (received_camera_measurement_)
        {
            model_->LOSFromGimbalAndCamera(x_(STATE_PITCH), x_(STATE_YAW), y_);
        }
        
    }
    
    received_gimbal_measurement_ = true;
    
    ukf.setMeasurementModel(model::stateToMeasurementGimbal);
    
    ukf.setMeasurementNoise(
            R.block(MEASUREMENT_GIMBAL_PITCH, MEASUREMENT_GIMBAL_PITCH, 5, 5));
    Vector11d u;
    u.setZero();
    ukf.predict(u);
    ukf.update(y_.segment(MEASUREMENT_GIMBAL_PITCH, 5));
}

void updateLRFMeasurement(double range)
{
    y_(MEASUREMENT_RANGE) = range;
    
    // Initialize range state on first reading
    if (not received_lrf_measurement_)
    {
        VectorXd x;
        x=ukf.getState();
        x(STATE_RANGE)=range;
        x(STATE_RANGE_RATE)=0;
        ukf.setState(x);
        ukf.setProcessFlags(RECEIVED_LRF, true);
//
//    P_(STATE_RANGE, STATE_RANGE) = initial_covariance_(STATE_RANGE);
//    P_(STATE_RANGE_RATE, STATE_RANGE_RATE) = initial_covariance_(STATE_RANGE_RATE);
        model_->initializeRangeAndRangeRate(range);
    }
    
    received_lrf_measurement_ = true;
    
    ukf.setMeasurementModel(model::stateToMeasurementLRF);
    ukf.setMeasurementNoise(R.block(MEASUREMENT_RANGE, MEASUREMENT_RANGE, 1, 1));
    Vector11d u;
    u.setZero();
    ukf.predict(u);
    ukf.update(y_.segment(MEASUREMENT_RANGE, 1));
}


void initializeUKF()
{
    // Initial covariance
    std::string path = getenv("COLCON_PREFIX_PATH");
    path.append("/eyeit/share/eyeit/params/los_estimator.params.yaml");
    YAML::Node paramNode = YAML::LoadFile(path);
    
    Vector11d P;
    P(STATE_RANGE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["range"].as<double>();
    P(STATE_RANGE_RATE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["range_rate"].as<double>();
    P(STATE_PITCH) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["angle"].as<double>();
    P(STATE_PITCH_RATE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["angle_rate"].as<double>();
    P(STATE_YAW) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["angle"].as<double>();
    P(STATE_YAW_RATE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["angle_rate"].as<double>();
    P(STATE_GIMBAL_PITCH) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["gimbal_angle"].as<double>();
    P(STATE_GIMBAL_YAW) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["gimbal_angle"].as<double>();
    P(STATE_GYRO_X) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["gyro"].as<double>();
    P(STATE_GYRO_Y) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["gyro"].as<double>();
    P(STATE_GYRO_Z) = paramNode["/**"]["ros__parameters"]["los_estimator"]["initial_covariance"]["gyro"].as<double>();
    
    // Process noise
    Vector11d Q;
    Q(STATE_RANGE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["range"].as<double>();
    Q(STATE_RANGE_RATE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["range_rate"].as<double>();
    Q(STATE_PITCH) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["angle"].as<double>();
    Q(STATE_PITCH_RATE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["angle_rate"].as<double>();
    Q(STATE_YAW) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["angle"].as<double>();
    Q(STATE_YAW_RATE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["angle_rate"].as<double>();
    Q(STATE_GIMBAL_PITCH) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["gimbal_angle"].as<double>();
    Q(STATE_GIMBAL_YAW) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["gimbal_angle"].as<double>();
    Q(STATE_GYRO_X) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["gyro"].as<double>();
    Q(STATE_GYRO_Y) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["gyro"].as<double>();
    Q(STATE_GYRO_Z) = paramNode["/**"]["ros__parameters"]["los_estimator"]["process_noise"]["gyro"].as<double>();
    
    // Measurement noise
    R.setZero();
    R(MEASUREMENT_RANGE,
      MEASUREMENT_RANGE) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["range"].as<double>();
    R(MEASUREMENT_X1C,
      MEASUREMENT_X1C) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["vision"].as<double>();
    R(MEASUREMENT_Y1C,
      MEASUREMENT_Y1C) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["vision"].as<double>();
    R(MEASUREMENT_GIMBAL_PITCH,
      MEASUREMENT_GIMBAL_PITCH) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["gimbal_angle"].as<double>();
    R(MEASUREMENT_GIMBAL_YAW,
      MEASUREMENT_GIMBAL_YAW) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["gimbal_angle"].as<double>();
    R(MEASUREMENT_GYRO_X,
      MEASUREMENT_GYRO_X) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["gyro"].as<double>();
    R(MEASUREMENT_GYRO_Y,
      MEASUREMENT_GYRO_Y) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["gyro"].as<double>();
    R(MEASUREMENT_GYRO_Z,
      MEASUREMENT_GYRO_Z) = paramNode["/**"]["ros__parameters"]["los_estimator"]["measurement_noise"]["gyro"].as<double>();
    
    Vector11d X0;
    X0.setZero();
    VectorXb p_flags(2);
    p_flags.setZero();
    VectorXb m_flags(0);
    ukf.initialize(X0, P.asDiagonal(), Q.asDiagonal(), R, model::nonLinModel,p_flags,m_flags);
    
    max_range_rate_validity = paramNode["/**"]["ros__parameters"]["los_estimator"]["range_rate_validity"]["max_covariance"].as<double>();
    min_range_rate_validity = paramNode["/**"]["ros__parameters"]["los_estimator"]["range_rate_validity"]["min_covariance"].as<double>();
    max_range_validity = paramNode["/**"]["ros__parameters"]["los_estimator"]["range_validity"]["max_covariance"].as<double>();
    min_range_validity = paramNode["/**"]["ros__parameters"]["los_estimator"]["range_validity"]["min_covariance"].as<double>();
    
    
    //    x_(STATE_RANGE) = range;
//    x_(STATE_RANGE_RATE) = 0;

//    x_(STATE_PITCH_RATE) = 0;
//    x_(STATE_YAW_RATE) = 0;
//    x_(STATE_YAW) = los_yaw;
//    x_(STATE_PITCH) = los_pitch;
//
//    x_(STATE_GIMBAL_PITCH) = y_(MEASUREMENT_GIMBAL_PITCH);
//    x_(STATE_GIMBAL_YAW) = y_(MEASUREMENT_GIMBAL_YAW);
//    x_(STATE_GYRO_X) = y_(MEASUREMENT_GYRO_X);
//    x_(STATE_GYRO_Y) = y_(MEASUREMENT_GYRO_Y);
//    x_(STATE_GYRO_Z) = y_(MEASUREMENT_GYRO_Z);
}


void drawPlots()
{
    plt::plot(time_for_plot_range, range_measurements, {
            {"label",           "LRF Measurements"},
            {"color",           "blue"},
            {"marker",          "s"},
            {"markerfacecolor", "none"},
            {"markersize",      "1"},
            {"linewidth",       "1.5"},
            {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, estimated_range,
              {{"label",           "Estimated Range"},
               {"color",           "red"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, estimated_range_rate,
              {{"label",           "Estimated Range Rate"},
               {"color",           "green"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot_gps_los, gps_los_range,
              {{"label",     "gps_los Range"},
               {"color",     "magenta"},
               {"linewidth", "1.0"},
               {"linestyle", "--"},});
    
    plt::plot(time_for_plot_gps_los, gps_los_range_rate,
              {{"label",     "gps_los Range rate"},
               {"color",     "black"},
               {"linewidth", "1.0"},
               {"linestyle", "--"},});
    
    plt::grid(true);
    //  plt::axis("equal");
    plt::title("Range and Range Rate");
    plt::draw();
    plt::legend();
    plt::xlabel("Time [sec]");
    plt::ylabel("Range [m]");
    plt::ylim(-50, 500);
    
    plt::figure(2);
    
    plt::plot(time_for_plot_gimbal, yaw_measurements, {
            {"label",           "Gimbal Yaw Measurements"},
            {"color",           "blue"},
            {"marker",          "s"},
            {"markerfacecolor", "none"},
            {"markersize",      "1"},
            {"linewidth",       "1.5"},
            {"linestyle",       "-"},});
    
    plt::plot(time_for_plot_gimbal, pitch_measurements,
              {{"label",           "Gimbal Pitch Measurements"},
               {"color",           "red"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, estimated_pitch,
              {{"label",           "Estimated Pitch"},
               {"color",           "green"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, estimated_yaw,
              {{"label",           "Estimated Yaw"},
               {"color",           "magenta"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    
    plt::plot(time_for_plot, estimated_pitch_rate,
              {{"label",           "Estimated Pitch Rate"},
               {"color",           "black"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, estimated_yaw_rate,
              {{"label",           "Estimated Yaw Rate"},
               {"color",           "cyan"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    
    plt::plot(time_for_plot_gimbal, gyro_y_measurements,
              {{"label",     "Gyro Y Measurements"},
               {"color",     "purple"},
               {"linewidth", "0.5"},
               {"linestyle", "-"},});
    
    
    plt::plot(time_for_plot_gimbal, gyro_z_measurements,
              {{"label",     "Gyro Z Measurements"},
               {"color",     "chocolate"},
               {"linewidth", "0.5"},
               {"linestyle", "-"},});
    
    plt::plot(time_for_plot_vision, center_x_measurements,
              {{"label",     "center x measurements"},
               {"color",     "orange"},
               {"linewidth", "0.5"},
               {"linestyle", "-"},});
    
    
    plt::plot(time_for_plot_vision, center_y_measurements,
              {{"label",     "center y measurements"},
               {"color",     "brown"},
               {"linewidth", "0.5"},
               {"linestyle", "-"},});
    
    plt::plot(time_for_plot_gps_los, gps_los_pitch,
              {{"label",     "gps_los pitch"},
               {"color",     "green"},
               {"linewidth", "1.0"},
               {"linestyle", "--"},});
    
    plt::plot(time_for_plot_gps_los, gps_los_yaw,
              {{"label",     "gps_los yaw"},
               {"color",     "magenta"},
               {"linewidth", "1.0"},
               {"linestyle", "--"},});
    
    plt::plot(time_for_plot_gps_los, gps_los_pitch_rate,
              {{"label",     "gps_los pitch rate"},
               {"color",     "black"},
               {"linewidth", "1.0"},
               {"linestyle", "--"},});
    
    plt::plot(time_for_plot_gps_los, gps_los_yaw_rate,
              {{"label",     "gps_los yaw rate"},
               {"color",     "cyan"},
               {"linewidth", "1.0"},
               {"linestyle", "--"},});
    
    plt::grid(true);
    //  plt::axis("equal");
    plt::title("Angles and Angle Rates");
    plt::draw();
    plt::legend();
    plt::xlabel("Time [sec]");
    plt::ylabel("Angle [rad]");
    plt::ylim(-3, 3);
    
    
    plt::figure(3);
    
    plt::plot(time_for_plot, range_covariance, {
            {"label",           "Range Covariance"},
            {"color",           "blue"},
            {"marker",          "s"},
            {"markerfacecolor", "none"},
            {"markersize",      "1"},
            {"linewidth",       "1.5"},
            {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, range_rate_covariance,
              {{"label",           "Range Rate Covariance"},
               {"color",           "red"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, pitch_covariance,
              {{"label",           "Pitch Covariance"},
               {"color",           "green"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, pitch_rate_covariance,
              {{"label",           "Pitch Rate Covariance"},
               {"color",           "magenta"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    
    plt::plot(time_for_plot, yaw_covariance,
              {{"label",           "Yaw Covariance"},
               {"color",           "black"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::plot(time_for_plot, yaw_rate_covariance,
              {{"label",           "Yaw Rate Covariance"},
               {"color",           "cyan"},
               {"marker",          "s"},
               {"markerfacecolor", "none"},
               {"markersize",      "1"},
               {"linewidth",       "1.5"},
               {"linestyle",       "-"},});
    
    plt::grid(true);
    //  plt::axis("equal");
    plt::title("Covariances");
    plt::draw();
    plt::legend();
    plt::xlabel("Time [sec]");
    plt::ylabel("Covariance");
    plt::ylim(-100, 100);
    plt::show(true);
}

void checkRangeRateValidity()
{
    if (not range_rate_valid)
    {
        if (covariance(STATE_RANGE_RATE) < min_range_rate_validity)
        {
            ukf.setProcessFlags(RANGE_RATE_VALID,true);
    
            range_rate_valid = true;
            std::cout << "Range rate is valid!" << std::endl;
        }
    }
        
        // Range rate is valid
    else
    {
        if (covariance(STATE_RANGE_RATE) > max_range_rate_validity)
        {
            ukf.setProcessFlags(RANGE_RATE_VALID, false);
    
            range_rate_valid = false;
            std::cout << "Range rate is invalid!!!!" << std::endl;
        }
    }
}


void checkRangeValidity()
{
    if (not range_valid)
    {
        if (covariance(STATE_RANGE) < min_range_validity)
        {
            range_valid = true;
            std::cout << "Range is valid!" << std::endl;
        }
    }
        
        // Range rate is valid
    else
    {
        if (covariance(STATE_RANGE) > max_range_validity)
        {
            range_valid = false;
            std::cout << "Range is invalid!!!!" << std::endl;
        }
    }
}


bool isLrfMeasurementValid(double current_lrf_measurement, double current_lrf_time)
{
    if (current_lrf_measurement == -1)
    {
        return false;
    }
    
    if (abs((current_lrf_measurement - previous_lrf_measurement) / (current_lrf_time - previous_lrf_time)) > 50)
    {
        return false;
    }
    
    return true;
}


void readFile(std::vector<std::vector<double>> &target, const std::string &file_name)
{
    std::vector<double> row;
    std::string line, word;
    
    std::fstream file(file_name, std::ios::in);
    if (file.is_open())
    {
        getline(file, line);// title row
        while (getline(file, line))
        {
            row.clear();
            std::stringstream str(line);
            while (getline(str, word, ','))
            {
                row.push_back(std::stod(word));
            }
            target.push_back(row);
        }
    }
    else
    {
        std::cout << "Could not open the file\n";
    }
}


int main()
{
    initializeUKF();
    
    Vector3d gyro_measurement(0, 0, 0);
//      gimbal goes crazy
    std::string raw_data_file_name = "/home/daniel/eyeit_tests/3.11/bags/rosbag2_2022_11_03-11_17_47/2022_11_03-11_17_47_merged.csv";
    std::string gps_los_file_name = "/home/daniel/eyeit_tests/3.11/bags/rosbag2_2022_11_03-11_17_47/2022_11_03-11_17_47_gps_los.csv";
    /*   good interception
       std::string raw_data_file_name = "/home/daniel/eyeit_tests/3.11/bags/rosbag2_2022_11_03-13_25_34/2022_11_03-13_25_34_merged.csv" ;
       std::string gps_los_file_name = "/home/daniel/eyeit_tests/3.11/bags/rosbag2_2022_11_03-13_25_34/2022_11_03-13_25_34_gps_los.csv";
       */
    
    // move between pi -> -pi
//    std::string raw_data_file_name = "/home/daniel/eyeit_tests/13.11/bags/rosbag2_2022_11_13-14_05_52/2022_11_13-14_05_52_merged.csv";
//    std::string gps_los_file_name = "/home/daniel/eyeit_tests/13.11/bags/rosbag2_2022_11_13-14_05_52/2022_11_13-14_05_52_gps_los.csv";
    
    std::vector<std::vector<double>> raw_data;
    readFile(raw_data, raw_data_file_name);
    
    if (!gps_los_file_name.empty())
    {
        std::vector<std::vector<double>> gps_los_data;
        readFile(gps_los_data, gps_los_file_name);
        for (auto &row_item: gps_los_data)
        {
            time_for_plot_gps_los.push_back(row_item[0]);
            gps_los_range.push_back(row_item[1]);
            gps_los_pitch.push_back(row_item[2]);
            gps_los_yaw.push_back(row_item[3]);
            gps_los_range_rate.push_back(row_item[4]);
            gps_los_pitch_rate.push_back(row_item[5]);
            gps_los_yaw_rate.push_back(row_item[6]);
        }
    }
    
    
    for (auto &row_item: raw_data)
    {
        ukf.setDt(row_item[2]);
        time_for_plot.push_back(row_item[1]);
        covariance = ukf.getCovariance().diagonal();
        
        
        range_covariance.push_back(covariance(STATE_RANGE));
        range_rate_covariance.push_back(covariance(STATE_RANGE_RATE));
        pitch_covariance.push_back(covariance(STATE_PITCH));
        pitch_rate_covariance.push_back(covariance(STATE_PITCH_RATE));
        yaw_covariance.push_back(covariance(STATE_YAW));
        yaw_rate_covariance.push_back(covariance(STATE_YAW_RATE));
        switch ((int) row_item[0])
        {
            case 1: // gimbal
                time_for_plot_gimbal.push_back(row_item[1]);
                pitch_measurement = (double) row_item[3];
                yaw_measurement = (double) helpers::azimuth2yaw(row_item[4]);
                pitch_measurements.push_back(pitch_measurement);
                yaw_measurements.push_back(yaw_measurement);
                gyro_measurement.x() = (double) row_item[5];
                gyro_measurement.y() = (double) row_item[6];
                gyro_measurement.z() = (double) row_item[7];
                updateGimbalMeasurement(pitch_measurement, yaw_measurement, gyro_measurement);
                
                gyro_y_measurements.push_back(gyro_measurement.y());
                gyro_z_measurements.push_back(gyro_measurement.z());
                break;
            
            case 2: // vision
                time_for_plot_vision.push_back(row_item[1]);
                center_x = (double) row_item[8];
                center_y = (double) row_item[9];
                center_x_measurements.push_back(center_x);
                center_y_measurements.push_back(center_y);
                updateVisionMeasurement(center_x, center_y);
                break;
            
            case 3: // lrf
                
                range_measurement = (double) row_item[10];
                double current_lrf_time = row_item[1];
                
                if (not isLrfMeasurementValid(range_measurement, current_lrf_time))
                {
                    break;
                }
                
                
                time_for_plot_range.push_back(current_lrf_time);
                range_measurements.push_back(range_measurement);
                
                range_measurement *= sin(state(STATE_PITCH) + center_y) /
                                     sin(state(STATE_PITCH)); //geometrical compensation of range(based on vision error)
                
                
                updateLRFMeasurement(range_measurement);
                
                previous_lrf_measurement = range_measurement;
                previous_lrf_time = current_lrf_time;
                break;
        }
        
        state = ukf.getState();
        
        
        estimated_range.push_back(state((STATE_RANGE)));
        estimated_range_rate.push_back(state((STATE_RANGE_RATE)));
        estimated_pitch.push_back(state((STATE_PITCH)));
        estimated_yaw.push_back(state((STATE_YAW)));
        estimated_pitch_rate.push_back(state((STATE_PITCH_RATE)));
        estimated_yaw_rate.push_back(state((STATE_YAW_RATE)));
        checkRangeValidity();
        checkRangeRateValidity();
        
        // TODO: Handle this properly
        if (range_measurement < 20 and range_measurement > 0 and estimated_range_rate.back() < 0)
        {
            break;
        }
        
    }
    
    drawPlots();
}

