#ifndef FUNCTIONS_ST_H
#define FUNCTIONS_ST_H

#include <vector>

// Linspace
std::vector<double> linspace_st(double start, double end, int num);

// Higher characteristic velocity (|u| + a)
double maior_velocidade_st(const std::vector<std::vector<double>>& U,
                           const std::vector<double>& a);

// 3x3 matrix Product
std::vector<std::vector<double>> produto_matriz_st(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B);

// Matrix-vector product: 3x3 * 3x1
std::vector<double> produto_matriz_vetor_st(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& v);

// Residual
double calcular_residuos_st(const std::vector<std::vector<double>>& Res,
                            double& res_1,
                            double& res_2,
                            double& res_3);

#endif // FUNCTIONS_ST_H
