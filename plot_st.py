import numpy as np
import matplotlib.pyplot as plt

# Function to read files
def ler_dados(nome_arquivo):
    return np.loadtxt(nome_arquivo)

# Reading the output files
x_p,   p   = ler_dados("pressure.txt").T
x_rho, rho = ler_dados("density.txt").T
x_u,   u   = ler_dados("velocity.txt").T

# Plot
fontSize       = 14
fontSizeLegend = 12
darkBlue  = (0.0,    0.129,  0.4784)
darkRed   = (0.7176, 0.0705, 0.207)
darkGreen = (0.0,    0.392,  0.196)

plt.rc('font', family='serif')

fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Pressure
axs[0].plot(x_p, p, marker='none', markersize=3, linestyle='-', color=darkBlue)
axs[0].set_title("Pressure", fontsize=fontSize)
axs[0].set_xlabel("$x$", fontsize=fontSize)
axs[0].set_ylabel("$p$", fontsize=fontSize)
axs[0].grid(False)

# Density
axs[1].plot(x_rho, rho, marker='none', markersize=3, linestyle='-', color=darkRed)
axs[1].set_title("Density", fontsize=fontSize)
axs[1].set_xlabel("$x$", fontsize=fontSize)
axs[1].set_ylabel(r"$\rho$", fontsize=fontSize)
axs[1].grid(False)

# Velocity
axs[2].plot(x_u, u, marker='none', markersize=3, linestyle='-', color=darkGreen)
axs[2].set_title("Velocity", fontsize=fontSize)
axs[2].set_xlabel("$x$", fontsize=fontSize)
axs[2].set_ylabel("$u$", fontsize=fontSize)
axs[2].grid(False)

plt.suptitle("Shock Tube — Steger & Warming", fontsize=fontSize+2, fontweight='bold')
plt.tight_layout()
plt.savefig("painel_shock_tube.png", dpi=150)
plt.show()

# Individual Plots

# Pressure
plt.figure(figsize=(8, 5))
plt.plot(x_p, p, marker='none', markersize=3, linestyle='-', color=darkBlue)
plt.title("Pressure", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$p$", fontsize=fontSize)
plt.grid(False)
plt.tight_layout()
plt.savefig("pressao_st.png", dpi=150)
plt.show()

# Density
plt.figure(figsize=(8, 5))
plt.plot(x_rho, rho, marker='none', markersize=3, linestyle='-', color=darkRed)
plt.title("Density", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel(r"$\rho$", fontsize=fontSize)
plt.grid(False)
plt.tight_layout()
plt.savefig("densidade_st.png", dpi=150)
plt.show()

# Velocity
plt.figure(figsize=(8, 5))
plt.plot(x_u, u, marker='none', markersize=3, linestyle='-', color=darkGreen)
plt.title("Velocity", fontsize=fontSize)
plt.xlabel("$x$", fontsize=fontSize)
plt.ylabel("$u$", fontsize=fontSize)
plt.grid(False)
plt.tight_layout()
plt.savefig("velocidade_st.png", dpi=150)
plt.show()
