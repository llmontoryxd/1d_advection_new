#include "DG.h"
#include "algebra.h"

dg::dg(int _k, int _N, double _L, double _tf, double _dt, bool _CFL_use, double _a, double _CFL = 1) {
    k = _k;
    N = _N;
    L = _L;
    tf = _tf;
    dt = _dt;
    x = Eigen::VectorXd::Zero(N+1);
    xc = Eigen::VectorXd::Zero(N);
    dx = L / N;
    a = _a;
    CFL_use = _CFL_use;
    CFL = _CFL;
    if (CFL_use == true) {
        dt = CFL * dx / a;
    }

    for (int i = 0; i < N+1; i++) { 
        x(i) = dx * i;
    }

    for (int i = 0; i < N; i++) {
        xc(i) = dx * (i+0.5);
    }

    log.open("log.txt");

    if (log.is_open()) {
        std::cout << "Log file was created" << std::endl;
    } else {
        std::cout << "Error: log file was not created" << std::endl;
    }

    log << "k = " << k << ", N = " << N << ", L = " << L << ", tf = " << tf << ", dt = " << dt << ", a = " << a << ", CFL_use = " << CFL_use << ", CFL = " << CFL << std::endl;


    log << "Allocating memory for M, K, U, G..." << std::endl;
    M = Eigen::MatrixXd::Zero(k, k);
    K = Eigen::MatrixXd::Zero(k, k);
    U = Eigen::MatrixXd::Zero(N, k);
    G = Eigen::MatrixXd::Zero(N, k);
    u = Eigen::VectorXd::Zero(N);
    log << "Allocating memory completed!" << std::endl << std::endl;

}

dg::~dg() {
    log.close();
}

double dg::calc_analytical(const int i) {
    return u0(x(i) - a * tf);
}

double dg::g(const int i_left, const int i_right) {
    //if (a >= 0) return a*x(i_left);
    //else return a*x(i_right);

    if (a >= 0) return calc_analytical(i_left);
    else return calc_analytical(i_right);

    //if (a >= 0) return u(i_left);
    //else return u(i_right);
}

double dg::u0(const double x) {
    return sin((4*M_PI*x)/L);
}

void dg::initial_condition() {
    log << "Working with initial conditions..." << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            U(i, j) = get_initial_u(k, j);
        }
    }
    log << U << std::endl;
    log << "Working with initial conditions completed!" << std::endl << std::endl;
}

double dg::get_initial_u(const int n, const int j) {
    Eigen::ArrayXd roots(n);
    Eigen::ArrayXd weights(n);

    roots = get_roots(n);
    weights = get_weights(n, roots);

    double ua0 = 0;
    for (int i = 0; i < n; i++) {
        ua0 += weights(i) * u0(roots(i)) * Pn(j, roots(i)) * roots(i);
    }

    return ua0;
}

void dg::get_M() {
    log << "Computing M..." << std::endl;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            if (i == j) {
                M(i, j) = get_M_el(k, i, j);    
            }
        }
    }  
    log << M << std::endl;
    log << "Computing M completed!" << std::endl << std::endl;
}

void dg::get_K() {
    log << "Computing K..." << std::endl;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            K(i, j) = get_K_el(k, i, j);
        }
    }  
    log << K << std::endl;
    log << "Computing K completed!" << std::endl << std::endl;  
}

void dg::get_G() {
    //log << "Computing G..." << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            G(i, j) = g(i+1, i+2) + pow(-1, j) * g(i, i+1);
        }
    }
    //log << G << std::endl;
    //log << "Computing G completed!" << std::endl << std::endl;
}

void dg::RK4() {
    for (int i = 0; i < N; i++) {
        Eigen::VectorXd k1, k2, k3, k4;
        k1 = get_kn(i, U.row(i));
        k2 = get_kn(i, U.row(i) + dt*k1.transpose()/2);
        k3 = get_kn(i, U.row(i) + dt*k2.transpose()/2);
        k4 = get_kn(i, U.row(i) + dt*k3.transpose());

        Eigen::VectorXd U_new;
        //log << k1 << std::endl << k2 << std::endl << k3 << std::endl << k4 << std::endl;
        U_new = U.row(i) + dt/6 * (k1.transpose() + 2*k2.transpose() + 2*k3.transpose() + k4.transpose());
        U.row(i) = U_new;
    }
    get_G();

    u = Eigen::VectorXd::Zero(N);
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < k; j++) {
            sum += U(i, j) * Pn(j, xc(i));
        }
        u(i) = sum;
    }
}

Eigen::VectorXd dg::get_kn(const int i, const Eigen::VectorXd f) {
    double C = 2/dx;
    //log << (C * a * MK(M, K) * f - C * a * M.transpose() * (G.row(i)).transpose()).transpose() << std::endl;
    //log << f << std::endl;
    //log << (G.row(i)).transpose() << std::endl;
    return C * a * MK(M, K) * f + C * a * M.transpose() * (G.row(i)).transpose();
}

void dg::solver() {
    initial_condition();
    get_M();
    get_K();
    get_G();
    double t = 0;
    t_now = t;

    log << "Preparations completed!!! Transition to the main cycle." << std::endl << std::endl;
    log << "Calculation by the Runge-Kutta method, 4th order..." << std::endl;
    while (t <= tf) {
        log << "t = " << t << ", max coeff = " << (U.cwiseAbs2()).maxCoeff() << std::endl;
        //t += dt;
        t_now = t;
        RK4();
        name = std::to_string(t);
        write_to_dat();
        t += dt;
    }

    u = Eigen::VectorXd::Zero(N);
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < k; j++) {
            sum += U(i, j) * Pn(j, xc(i));
        }
        u(i) = sum;
    }
}

void dg::write_to_dat() {
    std::ofstream results;
    string nname = "out/" + name + ".dat";
    results.open(nname);

    for (int i = 0; i < N; i++) {
        results << xc(i) << " " << u(i) << std::endl;
    }

    results.close();
}

void dg::save_to_vtk(const string name) {

}


