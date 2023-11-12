#pragma once
#include <vector>
#include <string>
#include <iostream>

using VD = std::vector<double>;

void ReadInput(std::string path, int& n, double& H, double& T0, double& dT, double& rho, double& Re, bool& method);

void Output(std::string path, const VD& l, const VD& r);