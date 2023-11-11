#pragma once
#include "in_out_file.h"


double mu_T(double T);

void CreateMesh(VD& y_nodes, VD& y_cells, int n, double H, double& dy);

void InitializeFields(VD& T, VD& U, VD& mu, VD& y, double T0, double dT);