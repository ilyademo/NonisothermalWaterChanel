#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include "Solve_file.h"

int main()
{
    //Create fields
    int n;
    VD T, U, mu, y_nodes, y_cells;
    double T0, dT, rho, Re, C, H, dy;
    bool interp_method;

    //Reading Input Data
    ReadInput("Input.txt", n, H, T0, dT, rho, Re, interp_method);

    //Calc dp/dx
    C = 12. * Re * pow(mu_T(T0), 2) / (rho * pow(H, 3));
    std::cout << "GradP: " << C << "\n";

    //Create nodes and cells
    CreateMesh(y_nodes, y_cells, n, H, dy);
    
    //Fill in T, mu in cells
    InitializeFields(T, U, mu, y_nodes, T0, dT);

    //Solve Equation
    SolveEq(U, mu, T, dy, C, T0, dT);

    //Calculate Average Velocity
    CalcAverageVel(U, y_cells, dy, H, Re, rho);

    //Calculate Viscosity standart deviation
    CalcViscDeviation(mu);
    
    //Output y = y(U)
    Output("u_y_cells.plt", U, y_cells);
}