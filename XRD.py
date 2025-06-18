from pymatgen.core import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Poscar

# POSCARファイルを読み込む
poscar = Poscar.from_file("albiteglass-shock-crystall-2_125000")
structure = poscar.structure
# XRDパターンを計算（Cu-Kα: λ = 1.5406 Å）
xrd_calc = XRDCalculator(wavelength="CuKa")

# 回折パターンの計算
xrd_pattern = xrd_calc.get_pattern(structure)

# プロット
plt.figure(figsize=(8, 4))
plt.plot(xrd_pattern.x, xrd_pattern.y, label="XRD Pattern")
plt.xlabel("2θ (degrees)")
plt.ylabel("Intensity (a.u.)")
plt.title("Simulated XRD Pattern")
plt.grid(True)
plt.savefig("albiteglass-shock-crystall-2_125000.png")
plt.legend()
plt.tight_layout()
plt.show()

