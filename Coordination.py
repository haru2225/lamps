import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
from collections import defaultdict

# xyzファイル読み込み
u = mda.Universe("albite-glass2.xyz", format='XYZ')

positions = u.atoms.positions
types = u.atoms.names  # 原子種の名前（例: 'Na', 'Al', 'Si', 'O'）

#カットオフ距離はRDF（Radial Distribution Function）の第一ピークの位置を基準に決める
# 元素ペアごとのカットオフ距離（Å）を辞書で指定
cutoff_dict = {
    ('Na', 'O'): 3.1,
    ('Al', 'O'): 2.5,
    ('Si', 'O'): 2.2,
    # 同じ元素同士の距離も必要なら追加
    ('Na', 'Na'): 2.5,
    ('Al', 'Al'): 2.5,
    ('Si', 'Si'): 3.0,
    ('O', 'O'): 2.5,
    # 他のペアは最大値を使うなど
}

def get_cutoff(type1, type2):
    # 並び替えたタプルで探す（Na-O と O-Na を同じとするため）
    pair = tuple(sorted((type1, type2)))
    return cutoff_dict.get(pair, 3.0)  # デフォルト3.0Å

coord_numbers = []

for i, pos_i in enumerate(positions):
    type_i = types[i]
    dist = distances.distance_array(pos_i.reshape(1, 3), positions)[0]
    
    neighbors = 0
    for j, d in enumerate(dist):
        if i == j:
            continue
        type_j = types[j]
        cutoff = get_cutoff(type_i, type_j)
        if d < cutoff:
            neighbors += 1
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

