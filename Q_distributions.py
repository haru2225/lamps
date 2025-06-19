from collections import Counter

datafile = "Al-O_bonds"  # LAMMPSのdataファイル名

# === Bondsセクションの抽出 ===
reading_bonds = False
center_counter = Counter()

with open(datafile) as f:
    for line in f:
        line = line.strip()
        if line.startswith("Bonds"):
            reading_bonds = True
            next(f)  # スキップ：セクション見出しの次行（空行 or ヘッダ）
            continue
        elif line.startswith("Angles") or line.startswith("Velocities") or line.startswith("Masses") or line.startswith("Atoms"):
            # 次のセクションに到達したら終了
            reading_bonds = False

        if reading_bonds:
            if line == "" or line.startswith("#"):
                continue
            tokens = line.split()
            if len(tokens) < 4:
                continue
            center_id = int(tokens[2])  # 3列目：中心原子ID
            center_counter[center_id] += 1

# === Q4だけ数える ===
q4_count = sum(1 for count in center_counter.values() if count == 4)

print(f"Q4 tetrahedra (center atoms with 4 bonds): {q4_count} atoms")

total_tetrahedra = len(center_counter)
q4_ratio = q4_count / total_tetrahedra * 100

print(f"Total tetrahedra (with ≥1 O): {total_tetrahedra}")
print(f"Q4 tetrahedra: {q4_count} atoms")
print(f"Q4 ratio: {q4_ratio:.2f}%")

