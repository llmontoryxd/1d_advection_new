#include "basis.h"

double Pn(const int n, const double x) {
    if (n == 0) return 1;
    else if (n == 1) return x;
    else return ((2*n - 1) * x * Pn(n-1, x) - (n - 1) * Pn(n-2, x)) / n;
}

double dPn(const int n, const double x) {
    if (n == 0) return 0;
    else if (n == 1) return 1;
    else return (n/(1-x*x)) * (Pn(n-1, x) - x * Pn(n, x));
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

Eigen::ArrayXd get_roots(const int n) {
    Eigen::ArrayXd roots(n);
    double h = 0.01;
    double atol = 1e-12;
    int maxSteps = 100;
    int i = 0;

    for (double x = -1.0; x <= 1.0; x = x + h) {
        double root = search_root(n, x, x + h, atol, maxSteps);
        if (root != 999) {
            roots(i) = root;
            i++;
        }
    }
    return roots;
}

double calculate_weight(const int n, const int i, Eigen::ArrayXd roots) {
    double x_i = roots(i);

    return 2 / ((1 - x_i * x_i) * dPn(n, x_i) * dPn(n, x_i));
}

Eigen::ArrayXd get_weights(const int n, Eigen::ArrayXd roots) {
    Eigen::ArrayXd weights(n);

    for (int i = 0; i < n; i++) {
        weights(i) = calculate_weight(n, i, roots);
    }

    return weights;
}

double get_K_el(const int n, const int a, const int b) {
    Eigen::ArrayXd roots(n);
    Eigen::ArrayXd weights(n);

    roots = get_roots(n);
    weights = get_weights(n, roots);

    double K_el = 0;
    for (int i = 0; i < n; i++) {
        K_el += weights(i) * Pn(a, roots(i)) * dPn(b, roots(i));
    }

    return K_el;
}

double get_M_el(const int n, const int a, const int b) {
    Eigen::ArrayXd roots(n);
    Eigen::ArrayXd weights(n);

    roots = get_roots(n);
    weights = get_weights(n, roots);

    double M_el = 0;
    for (int i = 0; i < n; i++) {
        M_el += weights(i) * Pn(a, roots(i)) * Pn(b, roots(i));
    }

    return M_el;    
}

