# Summary

[はじめに](README.md)

# 運用規則

- [クリーンルーム実装の手順](process/clean-room.md)
- [記述の規則と出典区分](process/writing-rules.md)
- [ライセンス](process/license.md)

# 資料索引

- [公開資料の一覧](methods/lcao-openmx/sources/index.md)
- [公開資料の空白](methods/lcao-openmx/sources/gaps.md)

# 全体設計

- [記法・単位・用語](conventions/notation.md)
- [モジュール構成とデータの流れ](methods/lcao-openmx/design/architecture.md)
- [H(R)・S(R) 疎ブロック形式](methods/lcao-openmx/design/sparse-block-format.md)
- [並列化](methods/lcao-openmx/design/parallelization.md)

# モジュール

- [PAO 基底と実球面調和関数](methods/lcao-openmx/modules/pao-basis.md)
- [二中心積分：重なり・運動エネルギー](methods/lcao-openmx/modules/two-center-overlap-kinetic.md)
- [非局所擬ポテンシャル](methods/lcao-openmx/modules/nonlocal-pseudopotential.md)
- [中性原子ポテンシャルと射影子展開](methods/lcao-openmx/modules/neutral-atom.md)
- [実空間グリッド・Poisson・交換相関](methods/lcao-openmx/modules/grid-poisson-xc.md)
- [k 点と Bloch 和](methods/lcao-openmx/modules/kpoints-bloch.md)
- [一般化固有値問題と S の正則化](methods/lcao-openmx/modules/generalized-eigen.md)
- [密度行列と電荷ミキシング](methods/lcao-openmx/modules/density-mixing.md)
- [全エネルギーと力](methods/lcao-openmx/modules/energy-force.md)
- [非共線・スピン軌道計算のエネルギー・力・応力](methods/lcao-openmx/modules/noncollinear-energy-force.md)

# 発展的な機能

- [O(N) 法](methods/lcao-openmx/modules/on-methods.md)
- [DFT+U・制約 DFT・磁気的な機能](methods/lcao-openmx/modules/dftu-magnetism.md)
- [ESM 法・外部電場・電荷ドープ・仮想原子・内殻準位](methods/lcao-openmx/modules/esm-field-corehole.md)
- [厳密交換・ハイブリッド汎関数・分散力補正](methods/lcao-openmx/modules/exx-vdw.md)
- [NEGF による電子輸送](methods/lcao-openmx/modules/negf.md)
- [構造最適化・分子動力学・NEB](methods/lcao-openmx/modules/geometry-md.md)

# 解析機能

- [バンド・状態密度・電荷解析・交換相互作用・Berry 位相・光学応答](methods/lcao-openmx/modules/analysis.md)
- [Wannier 関数・バンドアンフォールディング・局在軌道解析](methods/lcao-openmx/modules/wannier-unfolding.md)

# テスト

- [テストの階層と fixture](testing/strategy.md)
