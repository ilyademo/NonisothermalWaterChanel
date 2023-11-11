#include <fstream>
#include <vector>
#include <sstream>

using VD = std::vector<double>;

void ReadInput(std::string path, int& n, double& H, double& T0, double& dT, double& rho, double& Re, bool& method) {
    std::ifstream inp(path);
    double tmp;
    VD input;
    std::string str{ "" };
    if (inp.is_open()) {
        do {
            if (!inp.eof()) {
                inp >> tmp;
                input.push_back(tmp);
            }
        } while (getline(inp, str));
    }
    inp.close();
    n = input[0];
    H = input[1];
    T0 = input[2];
    dT = input[3];
    rho = input[4];
    Re = input[5];
    method = input[6];
}

void Output(std::string path, const VD& l, const VD& r) {
    std::ofstream out(path);
    if (out.is_open()) {
        for (auto i = 0; i < l.size(); i++) {
            out << l[i] << " " << r[i] << "\n";
        }
    }
    out.close();
}