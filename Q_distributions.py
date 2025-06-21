from pymatgen.io.vasp import Poscar
from pymatgen.analysis.local_env import CutOffDictNN
from collections import defaultdict

# --- POSCAR読み込み ---
structure = Poscar.from_file("DATA_POSCAR").structure

# --- カットオフ距離の指定 ---
nn_finder = CutOffDictNN(cut_off_dict={
    ("Si", "O"): 2.0,
    ("Al", "O"): 2.0,
    ("O", "Si"): 2.0,
    ("O", "Al"): 2.0,
})

# --- サイトとインデックスのマッピング（逆引き辞書） ---
site_to_index = {site: i for i, site in enumerate(structure)}

# --- Qn分布を記録する辞書 ---
Q_counts_Si = defaultdict(int)
Q_counts_Al = defaultdict(int)

# --- 四面体ごとにQn数を計算 ---
for i, site in enumerate(structure):
    symbol = site.specie.symbol
    if symbol not in ["Si", "Al"]:
        continue

    # 近傍O原子の取得
    neighbors = nn_finder.get_nn_info(structure, i)
    O_neighbors = [n for n in neighbors if n["site"].specie.symbol == "O"]

    qn = 0
    for O_info in O_neighbors:
        O_site = O_info["site"]
        try:
            O_index = site_to_index[O_site]
        except KeyError:
            # O_siteが同じ座標でもインスタンス違いで一致しないときは最近接探索で解決
            O_index = structure.get_sites_in_sphere(O_site.coords, 0.01)[0][0].index

        # O周りの近傍を再探索（このOが何個のSi/Alと結合しているか）
        O_neighbors2 = nn_finder.get_nn_info(structure, O_index)
        bonded_SiAl_count = sum(1 for n2 in O_neighbors2 if n2["site"].specie.symbol in ["Si", "Al"])
        if bonded_SiAl_count >= 2:
            qn += 1  # 架橋酸素と判定

    if symbol == "Si":
        Q_counts_Si[qn] += 1
    elif symbol == "Al":
        Q_counts_Al[qn] += 1

# --- 出力 ---
print("SiのQn分布:")
for qn in sorted(Q_counts_Si.keys(), reverse=True):
    print(f"  Q{qn}: {Q_counts_Si[qn]} 個")

print("\nAlのQn分布:")
for qn in sorted(Q_counts_Al.keys(), reverse=True):
    print(f"  Q{qn}: {Q_counts_Al[qn]} 個")

