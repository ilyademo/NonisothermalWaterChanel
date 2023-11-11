#pragma once
#include "Initialize_file.h"

void SolveEq(VD& u, const VD& mu, const VD& T, double dy, double gradP, double T0, double dT, bool m = true);

void CalcAverageVel(const VD& U, const VD& y, double dy, double H, double Re, double rho);

void CalcViscDeviation(const VD& mu);