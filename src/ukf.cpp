#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
    // a flag
    is_initialized_ = false;

    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    n_x_ = 5;
    x_ = VectorXd(n_x_);

    // initial augmented state vector
    n_aug_ = 7;
    x_aug_ = VectorXd(n_aug_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    //create sigma point matrix
    n_sigma_ = 2*n_aug_ + 1;
    P_aug_ = MatrixXd(n_aug_, n_aug_);
    A_ = MatrixXd(n_aug_, n_aug_);

    // The Xsig_aug will be filled out on the fly
    Xsig_aug_ = MatrixXd(n_aug_, n_sigma_);

    // The Xsig_pred_ for keeping predictions for the sigma points
    Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2; // TODO: adjust this

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 2; // TODO: adjust this

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

    /**
    TODO:

    Complete the initialization. See ukf.h for other member properties.

    Hint: one or more values initialized above might be wildly off...
    */

    // Radar covariance matrix
    L_ = MatrixXd(2, 2);
    L_ << std_laspx_*std_laspx_, 0,
          0,                     std_laspy_*std_laspy_;

    R_ = MatrixXd(3, 3);
    R_ << std_radr_*std_radr_, 0,                       0,
          0,                   std_radphi_*std_radphi_, 0,
          0,                   0,                       std_radrd_*std_radrd_;

    weights_ = VectorXd(n_sigma_);
    weights_(0) = lambda_ / (lambda_+n_aug_);
    for (int i = 1; i < n_sigma_; i++)
        weights_(i) = 0.5 / (n_aug_+lambda_);

    //Init NIS values
    NIS_radar_= 0.0;
    NIS_lidar_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if(!is_initialized_) {
        //set the previous_timestamp_
        previous_timestamp_ = meas_package.timestamp_;

        if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho = meas_package.raw_measurements_(0);
            double phi = meas_package.raw_measurements_(1);
            // NOTE: meas_package.raw_measurements_(2) is rho_dot, but we don't need it here
            x_ << cos(phi) * rho,
                  sin(phi) * rho,
                  0,
                  0,
                  0;
        }
        else if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
            double px = meas_package.raw_measurements_(0);
            double py = meas_package.raw_measurements_(1);
            x_ << px,
                  py,
                  0,
                  0,
                  0;
        }
        else
            cout << "I don't recognize sensor type: " << meas_package.sensor_type_ << endl;

        // Done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    //This is where the magic happens
    Prediction(dt);

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
        UpdateRadar(meas_package);
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
        UpdateLidar(meas_package);
    else
        cout << "I don't recognize sensor type: " << meas_package.sensor_type_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
    /**
    TODO:

    Complete this function! Estimate the object's location. Modify the state
    vector, x_. Predict sigma points, the state, and the state covariance matrix.
    */

    //I like it when the method names signify what they're actually doing
    GenerateSigmaPoints();

    PredictSigmaPoints(dt);

    PredictMeanAndCovariance();
}

void UKF::GenerateSigmaPoints() {

    x_aug_.head(n_x_) = x_;
    x_aug_.tail(n_aug_ - n_x_) << 0, 0;

    //update the augmented covariance matrix
    P_aug_.topLeftCorner(n_x_, n_x_) = P_;
    P_aug_.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) << std_a_*std_a_, 0,
                                                          0,             std_yawdd_*std_yawdd_;
    //calculate square root of P_aug_
    A_ = P_aug_.llt().matrixL();

    //set first column of sigma point matrix
    Xsig_aug_.col(0) = x_aug_;

    //set remaining sigma points
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * A_.col(i);
        Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * A_.col(i);
    }
}

void UKF::PredictSigmaPoints(double dt) {
    VectorXd col_handle;
    for (int i = 0; i < n_sigma_; i++) {
        col_handle = Xsig_aug_.col(i);
        Xsig_pred_.col(i) = PredictOne(col_handle, dt);
    }
}

VectorXd UKF::PredictOne(VectorXd& x_aug, double dt) {

    // NOTE: we don't really need positions to predict the sigma points, but FYI:
    // double px = x_aug(0);
    // double py = x_aug(1);
    double v = x_aug(2);
    double phi = x_aug(3);
    double phi_dot = x_aug(4);
    double nu_a = x_aug(5);
    double nu_phi = x_aug(6);

    double dt2_by_2 = dt * dt / 2;

    VectorXd change(n_x_);
    if(fabs(phi_dot) > 0.0001) {
        double v_by_phi_dot = v / phi_dot;
        double phi_plus = phi + phi_dot*dt;
        change << v_by_phi_dot * (sin(phi_plus) - sin(phi)),
                  v_by_phi_dot * (-cos(phi_plus) + cos(phi)),
                  0,
                  phi_dot * dt,
                  0;
    } else {
        change << v * cos(phi) * dt,
                  v * sin(phi) * dt,
                  0,
                  0,
                  0;
    }

    VectorXd noise(n_x_);
    noise << dt2_by_2 * cos(phi) * nu_a,
             dt2_by_2 * sin(phi) * nu_a,
             dt * nu_a,
             dt2_by_2 * nu_phi,
             dt * nu_phi;

    return x_aug.head(n_x_) + change + noise;
}

void UKF::PredictMeanAndCovariance() {
    //predict state mean
    x_.fill(0.0);
    for(int i = 0; i < n_sigma_; i++)
        x_ += weights_(i) * Xsig_pred_.col(i);

    //predict state covariance matrix
    P_.fill(0.0);
    for(int i = 0; i < n_sigma_; i++) {
        VectorXd diff = (Xsig_pred_.col(i) - x_);
        diff(3) = NormalizeAngle(diff(3));

        P_ += weights_(i) * diff * diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
    z_ = meas_package.raw_measurements_;
    int n_z = 2;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

    //transform sigma points into measurement space
    for (int i = 0; i < n_sigma_; i++) {
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);

        Zsig(0, i) = p_x;
        Zsig(1, i) = p_y;
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < n_sigma_; i++)
        z_pred += weights_(i) * Zsig.col(i);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd diff = (Zsig.col(i) - z_pred);

        S += weights_(i) * diff * diff.transpose();
    }

    //add measurement noise covariance matrix
    S = S + L_;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        VectorXd z_diff = (Zsig.col(i) - z_pred);
        VectorXd x_diff = (Xsig_pred_.col(i) - x_);

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain
    MatrixXd K = Tc * S.inverse();

    //update state mean and covariance matrix
    VectorXd z_diff = z_ - z_pred;
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
    cout << x_ << endl;

    //NIS
    NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */

    z_ = meas_package.raw_measurements_;
    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

    //transform sigma points into measurement space
    for (int i = 0; i < n_sigma_; i++) {
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v   = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double rho = sqrt(p_x*p_x + p_y*p_y);
        Zsig(0, i) = rho;
        Zsig(1, i) = atan2(p_y, p_x);
        Zsig(2, i) = v*(p_x*cos(yaw) + p_y*sin(yaw)) / rho;
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < n_sigma_; i++)
        z_pred = z_pred + weights_(i) * Zsig.col(i);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = NormalizeAngle(z_diff(1));

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    S = S + R_;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    /*****************************************************************************
    *  UKF Update for Radar
    ****************************************************************************/
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < n_sigma_; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = NormalizeAngle(z_diff(1));

        //state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = NormalizeAngle(x_diff(3));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z_ - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();

    //calculate NIS
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}

double UKF::NormalizeAngle(double angle) {
    while(angle > M_PI)
        angle -= 2.*M_PI;
    while(angle < -M_PI)
        angle += 2.*M_PI;
    return angle;
}
