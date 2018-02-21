#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

/*
 * Copied from the Project 1 repo:
 * https://github.com/MTDzi/CarND-Extended-Kalman-Filter-Project
 *
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    VectorXd diff;
    VectorXd squared_diff;

    // Check the validity of the following inputs:
    if(estimations.size() != ground_truth.size() || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    // Accumulate squared residuals
    for(int i=0; i<estimations.size(); ++i) {
        diff = estimations[i] - ground_truth[i];
        squared_diff = diff.array() * diff.array();
        rmse += squared_diff;
    }

    // Calculate the mean
    rmse /= estimations.size();

    // Calculate the squared root
    rmse = rmse.array().sqrt();

    // Return the result
    return rmse;
}