#include "functions_st.h"

#include <cmath>
#include <algorithm>

// Linspace
std::vector<double> linspace_st(double start, double end, int num) {
    std::vector<double> v(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i)
        v[i] = start + i * step;
    return v;
}

// Higher characteristic velocity (|u| + a)
double maior_velocidade_st(const std::vector<std::vector<double>>& U,
                           const std::vector<double>& a) {
    double vel_max = 0.0;
    size_t n = U[0].size();
    for (size_t i = 0; i < n; ++i) {
        double ui = U[1][i] / U[0][i];
        double val = std::abs(ui) + a[i];
        vel_max = std::max(vel_max, val);
    }
    return vel_max;
}

// 3x3 matrix Product
std::vector<std::vector<double>> produto_matriz_st(
            const std::vector<std::vector<double>>& A,
            const std::vector<std::vector<double>>& B) {

    std::vector<std::vector<double>> C(3, std::vector<double>(3, 0.0));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// Matrix-vector product: 3x3 * 3x1
std::vector<double> produto_matriz_vetor_st(
            const std::vector<std::vector<double>>& A,
            const std::vector<double>& v) {

    std::vector<double> res(3, 0.0);
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            res[i] += A[i][k] * v[k];
    return res;
}

// Residual
double calcular_residuos_st(const std::vector<std::vector<double>>& Res,
                            double& res_1,
                            double& res_2,
                            double& res_3) {

    res_1 = 0.0; res_2 = 0.0; res_3 = 0.0;
    int N = static_cast<int>(Res[0].size());

    for (int i = 1; i < N - 1; ++i) {
        res_1 = std::max(res_1, std::abs(Res[0][i]));
        res_2 = std::max(res_2, std::abs(Res[1][i]));
        res_3 = std::max(res_3, std::abs(Res[2][i]));
    }
    return std::max({res_1, res_2, res_3});
}
