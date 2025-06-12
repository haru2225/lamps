import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
from collections import defaultdict

# xyzファイル読み込み
u = mda.Universe("albite-test.xyz", format='XYZ')

cutoff = 3.0  # カットオフ距離[Å]
positions = u.atoms.positions
types = u.atoms.names  # 原子種の名前（例: 'Na', 'Al', 'Si', 'O'）

coord_numbers = []
for i, pos_i in enumerate(positions):
    dist = distances.distance_array(pos_i.reshape(1, 3), positions)[0]
    neighbors = np.sum((dist > 0) & (dist < cutoff))
    coord_numbers.append(neighbors)

# ファイルに書き込み
with open("coordination.txt", "w") as f:
    # 各原子の配位数をファイルに出力
    for i, cn in enumerate(coord_numbers):
        line = f"Atom {i} ({types[i]}) coordination number: {cn}\n"
        print(line, end='')  # 画面にも出力
        f.write(line)
    
    # 原子種ごとの平均配位数を計算してファイルに出力
    coord_sum = defaultdict(int)
    count = defaultdict(int)

    for t, cn in zip(types, coord_numbers):
        coord_sum[t] += cn
        count[t] += 1
    
    f.write("\nAverage coordination numbers per atom type:\n")
    print("\nAverage coordination numbers per atom type:")
    for t in coord_sum:
        avg_cn = coord_sum[t] / count[t]
        line = f"Type {t}: {avg_cn:.2f}\n"
        print(line, end='')
        f.write(line)

