#include "in_out_file.h"
#include "Initialize_file.h"
#include <fstream>

double mu_T(double T){
    return 3.3393e-11 * pow(T, 4) - 8.8614e-9 * pow(T, 3) + 9.3957e-7 * pow(T, 2) - 5.3115e-5 * T + 1.7579e-3;
}

void CreateMesh(VD& y_nodes, VD& y_cells, int n, double H, double& dy) {
    dy = H / (n - 1.);
    for (auto i = 0; i < n; i++) {
        y_nodes.push_back(dy * i);
        if (i != n - 1) {
            y_cells.push_back(dy / 2 + dy * i);
        }
    }
}

void InitializeFields(VD& T, VD& U, VD& mu, VD& y, double T0, double dT) {
    U.push_back(0);
    for (auto i = 0; i < y.size(); i++) {
        T.push_back(T0 - dT + 2 * dT * (y[i] - y[0]) / (y[y.size() - 1] - y[0]));
        U.push_back(0);
        if (i != 0) {
            T[i - 1] = (T[i] + T[i - 1]) / 2;
            mu.push_back(mu_T(T[i - 1]));
        }
    }
    T.erase(T.end() - 1);
    //Add fictive cells for interpolation
    T.insert(T.begin(), 2 * (T0 - dT) - T[0]);
    T.push_back(2 * (T0 + dT) - T[T.size() - 1]);
    mu.insert(mu.begin(), mu_T(T[0]));
    mu.push_back(mu_T(T[T.size() - 1]));
    Output("output_visc_cells.plt", T, mu);
}