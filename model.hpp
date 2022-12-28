//
// Created by daniel on 12/27/22.
//

#ifndef LOS_RATES_ESTIMATOR_MODEL_HPP
#define LOS_RATES_ESTIMATOR_MODEL_HPP

#include "UKF.hpp"

typedef Matrix<double, 23, 1> Vector23d;
typedef Matrix<double, 11, 1> Vector11d;
typedef Matrix<double, 8, 1> Vector8d;


enum State
{
    STATE_RANGE,
    STATE_RANGE_RATE,
    STATE_PITCH,
    STATE_PITCH_RATE,
    STATE_YAW,
    STATE_YAW_RATE,
    STATE_GIMBAL_PITCH,
    STATE_GIMBAL_YAW,
    STATE_GYRO_X,
    STATE_GYRO_Y,
    STATE_GYRO_Z,
};


enum Measurement
{
    MEASUREMENT_RANGE,
    MEASUREMENT_X1C,
    MEASUREMENT_Y1C,
    MEASUREMENT_GIMBAL_PITCH,
    MEASUREMENT_GIMBAL_YAW,
    MEASUREMENT_GYRO_X,
    MEASUREMENT_GYRO_Y,
    MEASUREMENT_GYRO_Z,
};

enum Sensor
{
    LRF,
    GIMBAL,
    CAMERA,
};

enum process_flags{
    RECEIVED_LRF,
    RANGE_RATE_VALID,
};


class model
{
public:
    model();
    
    ~model();
    
    static void stateToMeasurementVision(const MatrixXd &X, MatrixXd &Z, VectorXb flags);
    
    static void stateToMeasurementGimbal(const MatrixXd &X, MatrixXd &Z, VectorXb flags);
    
    static void stateToMeasurementLRF(const MatrixXd &X, MatrixXd &Z, VectorXb flags);
    
    static void nonLinModel(const MatrixXd &X, const VectorXd &u, MatrixXd &X_pred, double dt, VectorXb flags);
    
    static void LOSFromGimbalAndCamera(double &los_pitch, double &los_yaw, Vector8d y);
    
    
    static void initializeRangeAndRangeRate(double range);

private:
    static Matrix3d worldToCameraRotation(double gimbal_pitch, double gimbal_yaw);
    
    
    static void
    visionFromGimbalAndLos(double los_pitch, double los_yaw, double gimbal_pitch, double gimbal_yaw, double &x1c,
                           double &y1c);
    
};


#endif //LOS_RATES_ESTIMATOR_MODEL_HPP
