#include "ukf.h"
#include "Eigen/Dense"
#include "iostream"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  is_initialized_ = false;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  
  P_.fill(0.0);
  x_.fill(0.0);
  
  time_us_ = 0.0;

  weights_ = VectorXd(2*n_aug_+1);

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;

  for (int i=1; i<2*n_aug_+1; ++i) {  // 2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_){

    is_initialized_ = true;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER){

      //! Initialize state vector
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0,
            0,
            0;

      //! Initialize covariance matrix
      P_ << std_laspx_*std_laspy_,0,0,0,0,
            0, std_laspy_*std_laspy_,0,0,0,
            0,0,1,0,0,
            0,0,0,1,0,
            0,0,0,0,1;

    }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
 
      //! Initialize state vector
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double px = rho*cos(phi);
      double py = rho*sin(phi);
      x_ << px,
            py,
            0,
            0,
            0;

      //! Initialize covariance matrix
      P_ << std_radr_*std_radr_,0,0,0,0,
            0,std_radr_*std_radr_,0,0,0,
            0,0,std_radr_*std_radr_,0,0,
            0,0,0,std_radphi_,0,
            0,0,0,0,std_radphi_;
    }

  }else{
   
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000; //!convert to seconds
    Prediction(delta_t);
    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){    
       UpdateLidar(meas_package);
    }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
       UpdateRadar(meas_package);
    } 

  }
 
  time_us_ = meas_package.timestamp_;

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  //! Generate SIGMA points
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd p_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd XSig_aug = MatrixXd(n_aug_, 2*n_aug_+1); 

  //Augmented mean state
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //augmented covariance matrix
  p_aug.fill(0.0);
  p_aug.topLeftCorner(5,5) = P_;
  p_aug(5,5) = std_a_*std_a_;
  p_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = p_aug.llt().matrixL();

  //Create augmented sigma points
  XSig_aug.fill(0.0);
  XSig_aug.col(0) = x_aug;

  for(int i=0; i<n_aug_; ++i){
    XSig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    XSig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //predict sigma points
  for(int i=0; i<(2*n_aug_+1); ++i){
     // extract values for better readability
    double p_x = XSig_aug(0,i);
    double p_y = XSig_aug(1,i);
    double v = XSig_aug(2,i);
    double yaw = XSig_aug(3,i);
    double yawd = XSig_aug(4,i);
    double nu_a = XSig_aug(5,i);
    double nu_yawdd = XSig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  //predict state x
  x_.fill(0.0);
  for(int i=0; i<(2*n_aug_+1); ++i){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predict covariance matrix
  P_.fill(0.0);
  for(int i=0; i<(2*n_aug_+1); ++i){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
      
    while(x_diff(3) > M_PI) x_diff(3)-=2.*M_PI;
    while(x_diff(3) <-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  MatrixXd RLaser = MatrixXd (2,2);
  RLaser << (std_laspx_* std_laspx_), 0,
            0, (std_laspy_*std_laspy_);

  VectorXd ZCurr = VectorXd(2);
  ZCurr << meas_package.raw_measurements_[0],
           meas_package.raw_measurements_[1];

  MatrixXd Zsig = MatrixXd(2, 2*n_aug_+1);

  for(int i=0; i<(2*n_aug_+1); ++i){
   Zsig(0,i) = Xsig_pred_(0,i);
   Zsig(1,i) = Xsig_pred_(1,i);
  } 
  
  VectorXd z_pred = VectorXd(2);
  z_pred.fill(0.0);
  
  for(int i=0; i<(2*n_aug_+1); ++i){
    z_pred += weights_(i) * Zsig.col(i);
  }

  //! measurement covariance matrix
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);

  for(int i=0; i<(2*n_aug_+1); ++i){
    VectorXd Zdiff = Zsig.col(i) - z_pred;
    S += weights_(i) * Zdiff * Zdiff.transpose();
  }

  S = S + RLaser;
  
  //Crossscorrelation matrix
  MatrixXd Tc = MatrixXd(n_x_, 2);
  Tc.fill(0.0);

  for(int i=0; i<(2*n_aug_+1); ++i){
    //!residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //!state diference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();

  }
  
  //! Kalman gain
  MatrixXd K = Tc * S.inverse();
  
  //! residual
  VectorXd Z_diff = ZCurr - z_pred;

  //! update state
  x_ = x_ + K * Z_diff;

  //! update covariance matrix
  P_ = P_ - K * S * K.transpose();
    
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  
  int n_z = 3; //Radar dimensions
 
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;

  MatrixXd ZSig = MatrixXd(n_z, 2*n_aug_+1);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    ZSig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    ZSig(1,i) = atan2(p_y,p_x);                                // phi

    if(ZSig(0,i) < 1e-6){
       ZSig(0,i) = 1e-6;;
    }
    
    ZSig(2,i) = (p_x*v1 + p_y*v2) / ZSig(0,i); //sqrt(p_x*p_x + p_y*p_y);   // r_dot
    
  }

 VectorXd ZCurr = VectorXd(n_z);
 ZCurr << meas_package.raw_measurements_[0],
           meas_package.raw_measurements_[1],
           meas_package.raw_measurements_[2];
  
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  for(int i=0; i<(2*n_aug_+1); ++i){
    z_pred += weights_(i) * ZSig.col(i);
  }
  
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  
  for(int i=0; i<(2*n_aug_+1); ++i){
    VectorXd Zdiff = ZSig.col(i) - z_pred;

    // angle normalization
    while (Zdiff(1)> M_PI) Zdiff(1)-=2.*M_PI;
    while (Zdiff(1)<-M_PI) Zdiff(1)+=2.*M_PI;

    S += weights_(i) * Zdiff * Zdiff.transpose();

  }  

  S =  S + R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i=0; i<(2*n_aug_+1); ++i){
    //!residual
    VectorXd Zdiff = ZSig.col(i) - z_pred;
    
    // angle normalization
    while (Zdiff(1)> M_PI) Zdiff(1)-=2.*M_PI;
    while (Zdiff(1)<-M_PI) Zdiff(1)+=2.*M_PI;
    
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

  // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * Zdiff.transpose();

  }

  //! Kalman gain
  MatrixXd K = Tc * S.inverse();
  
  //! residual
  VectorXd Zdiff = ZCurr - z_pred;

  // angle normalization
  while (Zdiff(1)> M_PI) Zdiff(1)-=2.*M_PI;
  while (Zdiff(1)<-M_PI) Zdiff(1)+=2.*M_PI;

  // update state 
  x_ = x_ + K * Zdiff;  
  // update covariance
  P_ = P_ - K * S * K.transpose();
  
}
