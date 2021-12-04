#ifndef DG_H
#define DG_H

#include "basis.h"

class dg {
private:
// Variables
    // Mesh variables
    int k; // degree of the polynomial space
    int N; // number of partitions
    double L; // max coordinate value along the x-axis
    double tf; // final time
    double dt; // time step
    double dx; // coordinate step along the x-axis
    Eigen::VectorXd x; 
    Eigen::VectorXd xc;
    double t_now;

    // Equation variables
    double a; // transfer rate

    // DG method variables
    Eigen::MatrixXd M; // M matrix (k x k)
    Eigen::MatrixXd K; // K matrix (k x k)
    Eigen::MatrixXd U; // U matrix (N x k)
    Eigen::MatrixXd G; // G matrix (N x k)

    // Convegence variables
    bool CFL_use;
    double CFL;

    // Answer variables
    Eigen::VectorXd u; // what we need
    std::ofstream log;
    string name;


public:
// Methods
    // Constructor and destructor
    dg(int _k, int _N, double _L, double _tf, double _dt, bool _CFL_use, double _a, double _CFL);
    ~dg();

    // Flux function
    double g(const int, const int);

    // Initial condition
    void initial_condition();
    double u0(const double x);
    double get_initial_u(const int, const int);
    
    // Compute M, K, G
    void get_M();
    void get_K();
    void get_G();

    // Solver
    void RK4();
    Eigen::VectorXd get_kn(const int, const Eigen::VectorXd);
    void solver();

    // Writing results
    void write_to_dat();

    // Get analytical solution
    double calc_analytical(const int);

    // Save to vtk
    void save_to_vtk(const string);





};

#endif