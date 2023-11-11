#include "Initialize_file.h"

void SolveEq(VD& u, const VD& mu, const VD& T, double dy, double gradP, double T0, double dT, bool m = true) {
    size_t n = u.size();
    VD A(n), B(n), C(n), D(n), alph(n), bet(n);
    double mu_T0 = mu_T(T0);
    A[0] = 0;
    B[0] = 1;
    C[0] = 1;
    D[0] = 0;
    A[n - 1] = 1;
    B[n - 1] = 1;
    C[n - 1] = 0;
    D[n - 1] = 0;
    alph[0] = C[0] / B[0];
    bet[0] = D[0] / B[0];
    for (auto i = 1; i < n - 1; i++) {
        if (dT != 0) {
            if (m) {
                A[i] = (mu[i - 1] + mu[i]) / 2;
                C[i] = (mu[i + 1] + mu[i]) / 2;
            }
            else {
                A[i] = mu_T((T[i - 1] + T[i]) / 2);
                C[i] = mu_T((T[i + 1] + T[i]) / 2);
            }
            B[i] = -(C[i] + A[i]);
            D[i] = -gradP * pow(dy, 2);
        }
        else {
            A[i] = 1;
            B[i] = -2;
            C[i] = 1;
            D[i] = (-gradP * pow(dy, 2)) / mu_T0;
        }
        alph[i] = C[i] / (B[i] - A[i] * alph[i - 1]);
        bet[i] = (D[i] - A[i] * bet[i - 1]) / (B[i] - A[i] * alph[i - 1]);
    }
    bet[n - 1] = (D[n - 1] - A[n - 1] * bet[n - 2]) / (B[n - 1] - A[n - 1] * alph[n - 2]);
    u[n - 1] = bet[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        u[i] = bet[i] - alph[i] * u[i + 1];
    }
    //clear fictive U
    u.erase(u.end() - 1);
    u.erase(u.begin());
}

void CalcAverageVel(const VD& U, const VD& y, double dy, double H, double Re, double rho) {
    double res_num{ 0 }, res_analytic{ 0 };
    double tmp = mu_T(45);
    for (int i = 0; i < U.size() - 1; ++i) {
        res_num += (U[i] + U[i + 1]) / 2 * dy;
    }
    res_num /= H;
    res_analytic = Re * tmp / (rho * H);
    std::cout << "Numerical U: " << res_num << "\n";
    std::cout << "Analytical U: " << res_analytic << "\n";
    std::cout << "Deviation of U: " << abs((res_num - res_analytic) / res_analytic) * 100 << " %\n";
}

void CalcViscDeviation(const VD& mu) {
    double x_av{ 0 };
    double S{ 0 };
    for (auto i = 0; i < mu.size(); i++) {
        x_av += mu[i];
    }
    x_av /= mu.size() - 1;
    for (auto i = 0; i < mu.size(); i++) {
        S += pow((mu[i] - x_av), 2);
    }
    S /= mu.size() - 2;
    S = pow(S, 0.5);
    std::cout << "Standart deviation of visc = " << S << "\n";
}