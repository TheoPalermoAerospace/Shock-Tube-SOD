/*-----------------------------------------------------
--------------- SOD Shock Tube Problem ----------------
------------ Steger and Warming F-V Scheme ------------
-----------------------------------------------------*/

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include "functions_st.h"

int main() {

	// p, u and rho files
	std::ofstream pressure_data("pressure.txt");
	std::ofstream velocity_data("velocity.txt");
	std::ofstream density_data("density.txt");

	// Parameters
	const double gamma = 1.4;                         // Specific Heat Ratio
	const double CFL = 0.5;                           // CFL Condition
	const int x_i = 0;                                // First Node Coordinate
	const int x_f = 1;                                // Last Node Coordinate   
	const double L = 1.0;                             // Shock Tube Length
	const int N = 400;                                // Number of Nodes
	double delta_x = L/(N-1);                         // Mesh Spacing
	std::vector<double> x = linspace_st(x_i, x_f, N); // Mesh

	// Initializing Vectors
	std::vector<double> rho(N, 0.0);                                    // Density
	std::vector<double> u(N, 0.0);                                      // Velocity
	std::vector<double> a(N, 0.0);                                      // Sound Speed
	std::vector<double> e(N, 0.0);                                      // Total Energy
	std::vector<double> p(N, 0.0);                                      // Pressure
	std::vector<std::vector<double>> U(3, std::vector<double>(N, 0.0)); // State Vector

	// Steger & Warming FVS Matrices: A_pos, A_neg, F_pos, F_neg, lambda_pos, lambda_neg, P and inv(P)
	std::vector<std::vector<double>> A_pos(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> A_neg(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> F_pos(3, std::vector<double> (N, 0.0));
	std::vector<std::vector<double>> F_neg(3, std::vector<double> (N, 0.0));
	std::vector<std::vector<double>> lambda_pos(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> lambda_neg(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> P(3, std::vector<double> (3, 0.0));
	std::vector<std::vector<double>> P_inv(3, std::vector<double> (3, 0.0));

	// Initial Conditions
	for (int i = 0; i < N; ++i) {
		if (i <= N/2) {
			rho[i] = 1.0;
			p[i] = 1.0;
		} else {
			rho[i] = 0.125;
			p[i] = 0.1;
		}

		u[i] = 0.0;
		e[i] = p[i]/(gamma-1) + (0.5*rho[i]*pow(u[i], 2));
		a[i] = sqrt(gamma * p[i]/rho[i]);
	}

	// State Vector U
	for (int i = 0; i < N; ++i) {
		U[0][i] = rho[i];
		U[1][i] = rho[i] * u[i];
		U[2][i] = e[i];
	}

	// Auxiliar Vector
	std::vector<std::vector<double>> U_aux(3, std::vector<double>(N, 0.0));
	U_aux = U;

	// Time Marching Parameters
	double t_max = 0.2;
	double t = 0;

	std::cout << "\n----------------------------------------" << std::endl;
	std::cout << "-------- SOD Shock-Tube Problem --------" << std::endl;
	std::cout << "----------------------------------------" << std::endl;

	// Solver
	while (t < t_max) {

    // Calculate delta_t
    double vel_max = maior_velocidade_st(U, a);
    double delta_t = CFL*delta_x/vel_max;

    /*--------------------------------------
    ----- Flux Vector - Internal Nodes -----
    --------------------------------------*/
    for (int i = 0; i < N; ++i) {

        lambda_pos[0][0] = (u[i] + std::abs(u[i])) / 2.0;
        lambda_pos[1][1] = ((u[i]+a[i]) + std::abs(u[i]+a[i])) / 2.0;
        lambda_pos[2][2] = ((u[i]-a[i]) + std::abs(u[i]-a[i])) / 2.0;

        lambda_neg[0][0] = (u[i] - std::abs(u[i])) / 2.0;
        lambda_neg[1][1] = ((u[i]+a[i]) - std::abs(u[i]+a[i])) / 2.0;
        lambda_neg[2][2] = ((u[i]-a[i]) - std::abs(u[i]-a[i])) / 2.0;

        // P Matrix
        P[0][0] = 1.0;
        P[0][1] = U[0][i]/(2*a[i]);
        P[0][2] = -U[0][i]/(2*a[i]);
        P[1][0] = U[1][i]/U[0][i];
        P[1][1] = (U[0][i]/(2*a[i]))*((U[1][i]/U[0][i])+a[i]);
        P[1][2] = -(U[0][i]/(2*a[i]))*((U[1][i]/U[0][i])-a[i]);
        P[2][0] = pow((U[1][i]/U[0][i]), 2)/2;
        P[2][1] = (U[0][i]/(2*a[i]))*(0.5*pow((U[1][i]/U[0][i]), 2) + ((U[1][i]/U[0][i])*a[i]) + (pow(a[i], 2)/(gamma-1)));
        P[2][2] = -(U[0][i]/(2*a[i]))*(0.5*pow((U[1][i]/U[0][i]), 2) - ((U[1][i]/U[0][i])*a[i]) + (pow(a[i], 2)/(gamma-1)));

        // P_inv Matrix
        P_inv[0][0] = 1 - (0.5*(gamma-1)*pow((U[1][i]/U[0][i]), 2)/pow(a[i], 2));
        P_inv[0][1] = (gamma-1)*(U[1][i]/U[0][i])/pow(a[i], 2);
        P_inv[0][2] = -(gamma-1)/pow(a[i], 2);
        P_inv[1][0] = (1/(U[0][i]*a[i]))*(0.5*(gamma-1)*pow((U[1][i]/U[0][i]), 2) - ((U[1][i]/U[0][i])*a[i]));
        P_inv[1][1] = (a[i] - ((gamma-1)*(U[1][i]/U[0][i])))/(U[0][i]*a[i]);
        P_inv[1][2] = (gamma-1)/(U[0][i]*a[i]);
        P_inv[2][0] = -(1/(U[0][i]*a[i]))*(0.5*(gamma-1)*pow((U[1][i]/U[0][i]), 2) + ((U[1][i]/U[0][i])*a[i]));
        P_inv[2][1] = (a[i] + ((gamma-1)*(U[1][i]/U[0][i])))/(U[0][i]*a[i]);
        P_inv[2][2] = -(gamma-1)/(U[0][i]*a[i]);

        // A_pos = P * lambda_pos * P_inv
        std::vector<std::vector<double>> temp_pos = produto_matriz_st(P, lambda_pos);
        A_pos = produto_matriz_st(temp_pos, P_inv);

        // A_neg = P * lambda_neg * P_inv
        std::vector<std::vector<double>> temp_neg = produto_matriz_st(P, lambda_neg);
        A_neg = produto_matriz_st(temp_neg, P_inv);

        // F_pos[i] = A_pos * U[:,i]
        std::vector<double> Fpos_i = produto_matriz_vetor_st(A_pos, {U[0][i], U[1][i], U[2][i]});
        std::vector<double> Fneg_i = produto_matriz_vetor_st(A_neg, {U[0][i], U[1][i], U[2][i]});

        for (int j = 0; j < 3; ++j) {
            F_pos[j][i] = Fpos_i[j];
            F_neg[j][i] = Fneg_i[j];
        }
    }

    /*-------------------------------------------
    -------- Updating U - Internal Nodes --------
    -------------------------------------------*/
    for (int i = 1; i < N-1; ++i) {
        for (int j = 0; j < 3; ++j) {
            U_aux[j][i] = U[j][i] - (delta_t/delta_x) * ((F_pos[j][i] + F_neg[j][i+1]) - (F_pos[j][i-1] + F_neg[j][i]));
        }
    }

    /*------------------------------------------
    ------- Updating Primite Variables --------
    ------------------------------------------*/
    for (int i = 0; i < N; ++i) {
        rho[i] = U[0][i];
        u[i]   = U[1][i] / U[0][i];
        e[i]   = U[2][i];
        p[i]   = (e[i] - 0.5*rho[i]*pow(u[i], 2)) * (gamma-1);
        a[i]   = sqrt(gamma * p[i]/rho[i]);
    }

    std::cout << "Time Step: " << t << std::endl;
    t += delta_t;

    U = U_aux;

}

	// Writing p, u, rho
	for (int i = 0; i < N; ++i) {
    	pressure_data << x[i] << " " << p[i] << std::endl;
   		velocity_data << x[i] << " " << u[i]  << std::endl;
   		density_data << x[i] << " " << rho[i]  << std::endl;
	}

    std::cout << "\n---------------------------------------" << std::endl;
    std::cout << "---------- End of Simulation ----------" << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    pressure_data.close();
    velocity_data.close();
    density_data.close();
	
	return 0;
}
