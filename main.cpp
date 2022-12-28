#include <iostream>
#include "UKF.hpp"


void measurementModel(const MatrixXd &X, MatrixXd &Z){

}

int main()
{
    
    int dim_state = 7;
    int dim_measurements = 5;
    VectorXd x0(dim_state);// initial state
    VectorXd x(dim_state);// initial state
    x0.setZero();
    std::cout << "x0: \n" <<   x0 << std::endl;
    
    
    MatrixXd P0(dim_state, dim_state); // initial covariance
    MatrixXd P(dim_state, dim_state); // initial covariance
    P0 = MatrixXd::Identity(dim_state, dim_state)*10;
    std::cout << "P0: \n" <<   P0 << std::endl;
    
    
    MatrixXd Q(dim_state, dim_state); // Process noise
    Q.setZero();
    Q.diagonal() << 0.1,0.1,0.1,0.4,0.4,0.4,0.5;
    std::cout << "Q_: \n" <<   Q << std::endl;
    
    
    MatrixXd R(dim_measurements, dim_measurements); // measurements noise
    R.setZero();
    R.diagonal() << 0.1,0.2,0.3,0.4,0.5;
    std::cout << "R_: \n" <<   R << std::endl;
    
    
    void (*F)(const MatrixXd &X, const VectorXd &u,MatrixXd &X_pred); // process Model
    F =  &nonLinModel;
    std::cout << "F: \n" <<   F << std::endl;
    
    
    
    MatrixXd H(dim_state, dim_measurements); // measurement model
    H = MatrixXd(dim_measurements,dim_state ).setZero();
    H.diagonal() = VectorXd(dim_measurements).setOnes() * 1;
    std::cout << "H: \n" <<   H << std::endl;
    
    
    UKF ukf_(x0, P0, Q, R, &nonLinModel);
    
    VectorXd u(dim_state);
    u << 0, 0, 1, 0, 0, 0, 0;
    ukf_.predict(u);
    x=ukf_.getState();
    P=ukf_.getCovariance();
    
    VectorXd z(dim_measurements);
    z << 0,0,0.12,0,0;
    
    std::cout << "x_: \n" <<  x << std::endl;
    std::cout << "P_: \n" <<   P << std::endl;
    
    ukf_.update(z);
    
    x=ukf_.getState();
    P=ukf_.getCovariance();
    std::cout << "x_: \n" <<  x << std::endl;
    std::cout << "P_: \n" <<   P << std::endl;
    return 0;
}
