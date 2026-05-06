import numpy as np
import matplotlib.pyplot as plt


# Function to read files
def ler_dados(nome_arquivo):
    return np.loadtxt(nome_arquivo)

# Reading the output files (numerical solution)
x_p,   p_num   = ler_dados("pressure.txt").T
x_rho, rho_num = ler_dados("density.txt").T
x_u,   u_num   = ler_dados("velocity.txt").T

# Reading the output files (exact solution)
x_p_ex,   p_ex   = ler_dados("pressure_exact.txt").T
x_rho_ex, rho_ex = ler_dados("density_exact.txt").T
x_u_ex,   u_ex   = ler_dados("velocity_exact.txt").T

# Plot settings
fontSize       = 14
fontSizeLegend = 12
darkBlue  = (0.0, 0.129, 0.4784)
darkRed   = (0.7176, 0.0705, 0.207)
darkGreen = (0.0, 0.392, 0.196)
gray      = (0.4, 0.4, 0.4)

plt.rc('font', family='serif')

fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Pressure
axs[0].plot(x_p_ex, p_ex,   linestyle='--', color=gray,     lw=1.5, label='Exact Sol.')
axs[0].plot(x_p,    p_num,  linestyle='-',  color=darkBlue, lw=1.5, label='Numerical')
axs[0].set_title("Pressure", fontsize=fontSize)
axs[0].set_xlabel("$x$", fontsize=fontSize)
axs[0].set_ylabel("$p$", fontsize=fontSize)
axs[0].legend(fontsize=fontSizeLegend, frameon=False)
axs[0].grid(False)

# Density
axs[1].plot(x_rho_ex, rho_ex,  linestyle='--', color=gray,    lw=1.5, label='Exact Sol.')
axs[1].plot(x_rho,    rho_num, linestyle='-',  color=darkRed, lw=1.5, label='Numerical')
axs[1].set_title("Density", fontsize=fontSize)
axs[1].set_xlabel("$x$", fontsize=fontSize)
axs[1].set_ylabel(r"$\rho$", fontsize=fontSize)
axs[1].legend(fontsize=fontSizeLegend, frameon=False)
axs[1].grid(False)

# Velocity
axs[2].plot(x_u_ex, u_ex,   linestyle='--', color=gray,      lw=1.5, label='Exact Sol.')
axs[2].plot(x_u,    u_num,  linestyle='-',  color=darkGreen, lw=1.5, label='Numerical')
axs[2].set_title("Velocity", fontsize=fontSize)
axs[2].set_xlabel("$x$", fontsize=fontSize)
axs[2].set_ylabel("$u$", fontsize=fontSize)
axs[2].legend(loc='upper left', fontsize=fontSizeLegend, frameon=False)
axs[2].grid(False)

plt.tight_layout()
plt.savefig("painel_shock_tube.png", dpi=150)
plt.show()

# Individual Plots

# Pressure
plt.figure(figsize=(8, 5))
plt.plot(x_p_ex, p_ex,   linestyle='--', color=gray,     lw=1.5, label='Exact')
plt.plot(x_p,    p_num,  linestyle='-',  color=darkBlue, lw=1.5, label='Numerical')
plt.title("Pressure", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$p$", fontsize=fontSize)
plt.legend(fontsize=fontSizeLegend, frameon=False)
plt.grid(False)
plt.tight_layout()
plt.savefig("pressao_st.png", dpi=150)
plt.show()

# Density
plt.figure(figsize=(8, 5))
plt.plot(x_rho_ex, rho_ex,  linestyle='--', color=gray,    lw=1.5, label='Exact')
plt.plot(x_rho,    rho_num, linestyle='-',  color=darkRed, lw=1.5, label='Numerical')
plt.title("Density", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel(r"$\rho$", fontsize=fontSize)
plt.legend(fontsize=fontSizeLegend, frameon=False)
plt.grid(False)
plt.tight_layout()
plt.savefig("densidade_st.png", dpi=150)
plt.show()

# Velocity
plt.figure(figsize=(8, 5))
plt.plot(x_u_ex, u_ex,   linestyle='--', color=gray,      lw=1.5, label='Exact')
plt.plot(x_u,    u_num,  linestyle='-',  color=darkGreen, lw=1.5, label='Numerical')
plt.title("Velocity", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$u$", fontsize=fontSize)
plt.legend(loc='upper left', fontsize=fontSizeLegend, frameon=False)
plt.grid(False)
plt.tight_layout()
plt.savefig("velocidade_st.png", dpi=150)
plt.show()
