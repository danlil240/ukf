//
// Created by daniel on 12/26/22.
//

#ifndef LOS_RATES_ESTIMATOR_UKF_HPP
#define LOS_RATES_ESTIMATOR_UKF_HPP

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <memory>
#include <utility>

using namespace Eigen;
typedef Matrix<bool,Dynamic,1> VectorXb;

class UKF
{
public:
    UKF();
    
    UKF(VectorXd x0, MatrixXd P0, MatrixXd Q,
        MatrixXd R, void (*process_model_func)(const VectorXd &x, const VectorXd &u, VectorXd &x_pred, double dt, VectorXb flags),const VectorXb& process_flags,const VectorXb& measurements_flags);
    
    ~UKF();
    
    void initialize(VectorXd x0, MatrixXd P0, MatrixXd Q,
                    MatrixXd R,
                    void (*process_model_func)(const VectorXd &x, const VectorXd &u, VectorXd &x_pred, double dt, VectorXb flags),const VectorXb& process_flags,const VectorXb& measurements_flags);
    
    void predict(const VectorXd &u);
    
    void update(const VectorXd &z);
    
    void calcSigmaPoints(MatrixXd &X, VectorXd &W);
    
    VectorXd getState() const;
    
    MatrixXd getCovariance() const;
    
    void setProcessNoise(MatrixXd Q);
    
    void setMeasurementNoise(MatrixXd R);
    
    void setMeasurementModel(void (*func)(const MatrixXd &X, MatrixXd &Z,VectorXb flags));
    
    void setProcessFlags(int idx,bool value){
        process_flags_(idx)=value;
    }
    void setMeasurementsFlags(int idx,bool value){
        measurements_flags_(idx)=value;
    }
    
    void setState(VectorXd state){
        x_=std::move(state);
    }
    
    
    void setDt(double dt)
    { dt_ = dt; }

private:
    // Private member variables for the UKF class
    VectorXd x_;
    MatrixXd P_;
    MatrixXd Q_;
    MatrixXd R_;
    VectorXb process_flags_;
    VectorXb measurements_flags_;
    
    void (*processModel)(const VectorXd &x, const VectorXd &u, VectorXd &x_pred, double dt, VectorXb flags){};
    
    void (*measurementModel)(const MatrixXd &X, MatrixXd &Z,VectorXb flags){};
    
    double dt_{};
    
    
    void calcProcess(const MatrixXd &X, const VectorXd &u, MatrixXd &X_pred);
};

#endif //LOS_RATES_ESTIMATOR_UKF_HPP
