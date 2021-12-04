#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "basis.h"

Eigen::MatrixXd MK(Eigen::MatrixXd M, Eigen::MatrixXd K) {
    Eigen::MatrixXd Mt = M.transpose();
    return Mt*K;
}


#endif