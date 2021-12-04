#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

double Pn(const int n, const double x) {
    if (n == 0) return 1;
    else if (n == 1) return x;
    else return ((2*n - 1) * x * Pn(n-1, x) - (n - 1) * Pn(n-2, x)) / n;
}

double search_root(const int n, double a, double b, const double atol, const int maxSteps) { 
    double c;

    if (Pn(n, a) * Pn(n, b) <= 0) {
        int iter = 1;
        do {
            c = (a + b) / 2;
            if (Pn(n, a) * Pn(n, c) > 0) {
                a = c;
            } else if (Pn(n, a) * Pn(n, c) < 0) {
                b = c;
            } else if (Pn(n, c) == 0) {
                return c;
            }
            iter++;
        } while (fabs(a-b) >= atol && iter <= maxSteps);
        return c;
    } else {
        return 999;
    }
    
}

int main() {
    int n = 3;
    Eigen::ArrayXd roots(n);
    std::cout << roots << std::endl;
    double h = 0.01;
    double atol = 1e-12;
    int maxSteps = 10000;
    int i = 0;

    for (double x = -1.0; x <= 1.0; x = x + h) {
        double root = search_root(n, x, x + h, atol, maxSteps);
        if (root != 999) {
            //std::cout << root << std::endl; 
            roots(i) = root;
            i++;
        }
    }    
    std::cout << roots << std::endl;
}