//
// Created by daniel on 12/26/22.
//


#include "UKF.hpp"

#include <utility>

UKF::UKF()
{

};

UKF::UKF(VectorXd x0, MatrixXd P0, MatrixXd Q,
         MatrixXd R,
         void (*process_model_func)(const VectorXd &x, const VectorXd &u, VectorXd &x_pred, double dt, VectorXb flags),
         const VectorXb &process_flags, const VectorXb &measurements_flags)
        : x_(std::move(x0)), P_(std::move(P0)), Q_(std::move(Q)), R_(std::move(R))
{
    process_flags_ = process_flags;
    measurements_flags_ = measurements_flags;
    processModel = process_model_func;
    std::cout << "begin" << std::endl;
}


void UKF::initialize(VectorXd x0, MatrixXd P0, MatrixXd Q, MatrixXd R,
                     void (*process_model_func)(const VectorXd &x, const VectorXd &u, VectorXd &x_pred, double dt,
                                                VectorXb flags), const VectorXb &process_flags,
                     const VectorXb &measurements_flags)
{
    x_ = std::move(x0);
    P_ = std::move(P0);
    Q_ = std::move(Q);
    R_ = std::move(R);
    process_flags_ = process_flags;
    measurements_flags_ = measurements_flags;
    processModel = process_model_func;
}


UKF::~UKF()
{
    std::cout << "end" << std::endl;
}

void UKF::predict(const VectorXd &u)
{
    // Calculate sigma points
    MatrixXd X;
    VectorXd W;
    calcSigmaPoints(X, W);
    
    // Propagate sigma points through process model
    MatrixXd X_pred;
    calcProcess(X, u, X_pred);
    
    
    // Calculate mean and covariance of predicted state
    x_ = X_pred * W;
    
    P_ = X_pred * W.asDiagonal() * X_pred.transpose() + Q_;
}

void UKF::update(const VectorXd &z)
{
    // Calculate sigma points
    MatrixXd X;
    VectorXd W;
    calcSigmaPoints(X, W);
    
    // Propagate sigma points through measurement model
    MatrixXd Z;
    (*measurementModel)(X, Z, measurements_flags_);
    
    // Calculate mean and covariance of predicted measurement
    VectorXd z_pred = Z * W;
    
    // Calculate cross-covariance of state and measurement
    MatrixXd P_xz = (X.colwise() - x_) * (W.asDiagonal() * (Z.colwise() - z_pred).transpose());
    
    // Calculate Kalman gain
    MatrixXd K = P_xz * (Z * W.asDiagonal() * Z.transpose() + R_).inverse();
    
    // Update state and covariance
    x_ += K * (z - z_pred);
    P_ -= K * P_xz.transpose();
}

void UKF::calcSigmaPoints(Eigen::MatrixXd &X, Eigen::VectorXd &W)
{
    // Determine sigma point scaling factor
    double lambda = P_(0, 0) / (P_(0, 0) + R_(0, 0));
    
    // Calculate sigma points
    int n = x_.size();
    X = Eigen::MatrixXd::Zero(n, 2 * n + 1);
    X.col(0) = x_;
    X.block(0, 1, n, n) = (sqrt(n + lambda * R_(0, 0)) * Eigen::MatrixXd::Identity(n, n)).colwise() + x_;
    X.block(0, n + 1, n, n) = (-sqrt(n + lambda * R_(0, 0)) * Eigen::MatrixXd::Identity(n, n)).colwise() + x_;
    
    // Calculate weights
    W = Eigen::VectorXd::Ones(2 * n + 1);
    W(0) = lambda / (n + lambda * R_(0, 0));
    W.tail(2 * n).setConstant((1 - lambda) / (n + lambda * R_(0, 0)));
}

void UKF::calcProcess(const MatrixXd &X, const VectorXd &u, MatrixXd &X_pred)
{
    X_pred=MatrixXd::Zero(X.rows(),X.cols());
    for (int i = 0; i < X.cols(); i++)
    {
        VectorXd x_pred;
        processModel(X.col(i), u, x_pred, dt_, process_flags_);
        X_pred.col(i) = x_pred;
    }
}


VectorXd UKF::getState() const
{
    return x_;
}

MatrixXd UKF::getCovariance() const
{
    return P_;
    
}

void UKF::setProcessNoise(MatrixXd Q)
{
    Q_ = std::move(Q);
}

void UKF::setMeasurementNoise(MatrixXd R)
{
    R_.resize(R.rows(), R.cols());
    R_ = std::move(R);
    
}

void UKF::setMeasurementModel(void (*func)(const MatrixXd &X, MatrixXd &Z, VectorXb flags))
{
    measurementModel = func;
}
