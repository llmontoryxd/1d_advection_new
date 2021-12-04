#include "DG.h"

int main() {
    int k;
    int N;
    double L;
    double tf;
    double dt;
    bool CFL_use;
    double CFL;
    double a;

    ifstream fin;
    string line;
    string line1;

    fin.open("input_file.yaml");

    while (true) {
        getline(fin, line);
        if (line == "start input") break;
    }

    fin >> line1 >> k;
    fin >> line1 >> N;
    fin >> line1 >> L;
    fin >> line1 >> tf;
    fin >> line1 >> dt;
    fin >> line1 >> CFL_use;
    fin >> line1 >> a;
    fin >> line1 >> CFL;

    fin.close();

    dg solution(k, N, L, tf, dt, CFL_use, a, CFL);
    solution.solver();

}