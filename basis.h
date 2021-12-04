#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

using namespace std;


/* --------------------- Legendre Basis --------------------- */

double Pn(const int, const double); // P_n(x)
double dPn(const int, const double); // P'_n(x)
double search_root(const int, double, double, const double, const int); // Search root bisection method of P_n(x)
Eigen::ArrayXd get_roots(const int); // get all roots of P_n(x)
double calculate_weight(const int, const int); // get i-weight
Eigen::ArrayXd get_weights(const int, Eigen::ArrayXd); // get all weights
double get_K_el(const int, const int, const int); // get K_ab
double get_M_el(const int, const int, const int); // get M_a

#endif