//
// Created by daniel on 12/27/22.
//

#include "model.hpp"



model::model()
= default;

model::~model()
= default;



void model::nonLinModel(const MatrixXd &X, const VectorXd &u, MatrixXd &X_pred, double dt,VectorXb flags)
{
    X_pred=MatrixXd::Zero(11,23);
    Vector23d r = X.row(STATE_RANGE);
    Vector23d r_dot = X.row(STATE_RANGE_RATE);
    Vector23d pitch = X.row(STATE_PITCH);
    Vector23d pitch_rate = X.row(STATE_PITCH_RATE);
    Vector23d yaw = X.row(STATE_YAW);
    Vector23d yaw_rate = X.row(STATE_YAW_RATE);
    Vector23d gimbal_pitch = X.row(STATE_GIMBAL_PITCH);
    Vector23d gimbal_yaw = X.row(STATE_GIMBAL_YAW);
    Vector23d gyro_x = X.row(STATE_GYRO_X);
    Vector23d gyro_y = X.row(STATE_GYRO_Y);
    Vector23d gyro_z = X.row(STATE_GYRO_Z);
    
    // The acceleration of a vector in spherical coordinates:
    // https://en.wikipedia.org/wiki/Spherical_coordinate_system#Kinematics
    // Assuming no acceleration, each element (r, theta, phi) must be 0.
    // r_dot_dot, theta_dot_dot and phi_dot_dot can then be solved for
    
    Vector23d phi_dot = yaw_rate;
    Vector23d theta = pitch + Vector23d::Ones() * M_PI_2; // In the spherical kinematics equations, theta = 0 is up;
    Vector23d theta_dot = pitch_rate;
    
    Vector23d theta_dot_dot = Vector23d::Zero();
    Vector23d phi_dot_dot = Vector23d::Zero();
    Vector23d r_dot_dot = Vector23d::Zero();
    std::cout << "flags(RECEIVED_LRF): " <<   flags(RECEIVED_LRF) << std::endl;
    
    if (flags(RECEIVED_LRF))
    {
        for (int column = 0; column < 23; column++)
        {
            r_dot_dot(column) = r(column) * pow(theta_dot(column), 2) +
                                r(column) * pow(phi_dot(column), 2) * pow(sin(theta(column)), 2);
            
            if (flags(RANGE_RATE_VALID))//valid_range_rate_)
            {
                theta_dot_dot(column) = -2 * r_dot(column) * theta_dot(column) / r(column) +
                                        pow(phi_dot(column), 2) * sin(theta(column)) * cos(theta(column));
                phi_dot_dot(column) = -2 * r_dot(column) * phi_dot(column) / r(column) -
                                      2 * theta_dot(column) * phi_dot(column) * cos(theta(column)) /
                                      sin(theta(column));
            }
        }
        
        X_pred.row(STATE_RANGE) = r + r_dot * dt;
        X_pred.row(STATE_RANGE_RATE) = r_dot + r_dot_dot * dt;
    }
    
    X_pred.row(STATE_PITCH) = pitch + pitch_rate * dt;
    X_pred.row(STATE_PITCH_RATE) = pitch_rate + theta_dot_dot * dt;
    X_pred.row(STATE_YAW) = yaw + yaw_rate * dt;
    X_pred.row(STATE_YAW_RATE) = yaw_rate + phi_dot_dot * dt;
    X_pred.row(STATE_GIMBAL_PITCH) = gimbal_pitch + gyro_y * dt;
}


Matrix3d model::worldToCameraRotation(double gimbal_pitch, double gimbal_yaw)
{
    Matrix3d world_to_gimbal_y, world_to_gimbal_z;
    
    world_to_gimbal_z << cos(gimbal_yaw), -sin(gimbal_yaw), 0,
            sin(gimbal_yaw), cos(gimbal_yaw), 0,
            0, 0, 1;
    
    world_to_gimbal_y << cos(gimbal_pitch), 0, sin(gimbal_pitch),
            0, 1, 0,
            -sin(gimbal_pitch), 0, cos(gimbal_pitch);
    
    Matrix3d R_world_to_gimbal = world_to_gimbal_z * world_to_gimbal_y;
    
    Matrix3d gimbal_to_camera_y, gimbal_to_camera_z;
    
    gimbal_to_camera_y << 0, 0, 1, // Rotate 90 degrees pitch so that Z looks out from the camera
            0, 1, 0,
            -1, 0, 0;
    
    gimbal_to_camera_z << 0, 1, 0, // Rotate -90 degrees yaw so that X and Y are in the right directions
            -1, 0, 0,
            0, 0, 1;
    
    Matrix3d R_gimbal_to_camera = gimbal_to_camera_y * gimbal_to_camera_z;
    
    return R_world_to_gimbal * R_gimbal_to_camera;
}


void model::LOSFromGimbalAndCamera(double &los_pitch, double &los_yaw,Vector8d y)
{
    Vector3d target_position_in_image = {y(MEASUREMENT_X1C), y(MEASUREMENT_Y1C), 1};
    Vector3d los_direction_world_frame =
            worldToCameraRotation(y(MEASUREMENT_GIMBAL_PITCH), y(MEASUREMENT_GIMBAL_YAW)) * target_position_in_image;
    
    los_pitch = -atan2(los_direction_world_frame.z(), los_direction_world_frame.head<2>().norm());
    los_yaw = atan2(los_direction_world_frame.y(), los_direction_world_frame.x());
}


void model::visionFromGimbalAndLos(double los_pitch, double los_yaw, double gimbal_pitch, double gimbal_yaw,
                            double &x1c, double &y1c)
{
    // Calculate ENU LOS vector from angles
    Vector3d normalized_los_enu;
    double horizontal_los = cos(los_pitch);
    normalized_los_enu.x() = cos(los_yaw) * horizontal_los;
    normalized_los_enu.y() = sin(los_yaw) * horizontal_los;
    normalized_los_enu.z() = -sin(los_pitch);
    
    // LOS_direction = (R_world_to_gimbal * R_gimbal_to_camera) * vision_vector  ->  vision_vector = inv(R_world_to_camera) * LOS_direction
    Vector3d vision = worldToCameraRotation(gimbal_pitch, gimbal_yaw).transpose() * normalized_los_enu;
    
    x1c = vision.x();
    y1c = vision.y();
}


void model::initializeRangeAndRangeRate(double range)
{
    

//
//    P_(STATE_RANGE, STATE_RANGE) = initial_covariance_(STATE_RANGE);
//    P_(STATE_RANGE_RATE, STATE_RANGE_RATE) = initial_covariance_(STATE_RANGE_RATE);
}

void model::stateToMeasurementVision(const MatrixXd &X, MatrixXd &Z,VectorXb flags)
{
    Z = MatrixXd(2, 23);
    for (int i = 0; i < 23; i++)
    {
        double x1c, y1c;
        visionFromGimbalAndLos(X(STATE_PITCH, i), X(STATE_YAW, i), X(STATE_GIMBAL_PITCH, i), X(STATE_GIMBAL_YAW, i),
                               x1c, y1c);
        Z(0, i) = x1c;
        Z(1, i) = y1c;
    }
}


void model::stateToMeasurementGimbal(const MatrixXd &X, MatrixXd &Z,VectorXb flags)
{
    
    Z = MatrixXd(5, 23);
    
    Z.row(0) = X.row(STATE_GIMBAL_PITCH);
    Z.row(1) = X.row(STATE_GIMBAL_YAW);
    Z.row(2) = X.row(STATE_GYRO_X);
    Z.row(3) = X.row(STATE_GYRO_Y);
    Z.row(4) = X.row(STATE_GYRO_Z);
    
}


void model::stateToMeasurementLRF(const MatrixXd &X, MatrixXd &Z,VectorXb flags)
{
    Z = MatrixXd(1, 23);
    Z.row(0) = X.row(STATE_RANGE);
}

