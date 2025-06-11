import numpy as np
import matplotlib.pyplot as plt

# データ（密度と圧力）
density = np.array([1.3815, 1.5137, 1.6585, 1.8172, 1.9911, 2.1816])  # g/cm³
pressure = np.array([327629, 404494, 497842, 612622, 754935, 932757])  # atm

# 初期密度（最大体積時＝最も小さい密度と仮定）
rho0 = density[0]

# V/V0 = rho0 / rho
V_over_V0 = rho0 / density

# c^2 = dP/d(rho)
c_squared = np.gradient(pressure, density)

# プロット
fig, ax = plt.subplots(2, 1, figsize=(7, 10), sharex=True)

# 図A: Cold pressure vs V/V0
ax[0].plot(V_over_V0, pressure / 1e3, 'o-', color='blue')  # 単位：k atm
ax[0].set_ylabel('Pressure (k atm)')
ax[0].set_title('Cold Pressure vs V/V0')
ax[0].grid(True)

# 図B: Sound speed squared vs V/V0
ax[1].plot(V_over_V0, c_squared, 's-', color='red')
ax[1].set_xlabel('V/V0')
ax[1].set_ylabel('Sound Speed Squared (atm·cm³/g)')
ax[1].set_title('Sound Speed Squared vs V/V0')
ax[1].grid(True)

plt.tight_layout()
plt.show()

