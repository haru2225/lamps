# Albite Shock Simulation with LAMMPS

Albite結晶とAlbiteアモルファスに対する衝撃シュミレーションを行っています。
結晶・アモルファス構造の比較、XRD解析、Hugoniostat条件での挙動を検討しています。


・ ファイル一覧

- `albite-crystal-shock.lj`：結晶構造での衝撃シミュレーション
- `albite-amorphous-shock.lj`：アモルファス構造での衝撃シミュレーション 
- `albite-data31`, `albite-data40`：結晶とアモルファスの初期構造データ
- `albite-crystal-amorphous.lj`：結晶からアモルファス作成
- `albite-random-amorphous.lj`：ランダムな原子配置からアモルファス作成
- `albite-random-amorphous.xyz`：albite-random-amorphous.ljの最終的な構造データ
- `albite-crystal-amorphous.xyz`：albite-crystal-amorphous.ljの最終的な構造データ
- `XRD.py`：XRDパターン生成スクリプト（pymatgen使用）
- `XRD-albite-crystal.png`, `XRD-albite-amorphous.png`：XRD画像

・実装
以下のコマンドで実装

mpirun -np X ~/lmp_mpi -in  Filename

ex)mpirun -np 16 ~/lmp_mpi -in  albite-random-amorphous.lj


・その他:
結晶構造データはAmerican Mineralogist Crystal Structure Database から得ています。


