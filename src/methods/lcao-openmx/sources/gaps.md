# 公開資料の空白

公開資料（技術ノート・論文・マニュアル）だけでは再実装に足りない項目の一覧です。
各項目は、記述者が〔規約〕として補うか、このノートで〔設計〕として決め直します。

## 全体に関わる規約

| 項目 | 状況 | 対応 |
|---|---|---|
| 実球面調和関数の並び順・符号 | マニュアルは「実球面調和関数」とのみ記載 | [PAO 基底](../modules/pao-basis.md) に〔規約〕として記載済み |
| PAO・VPS ファイルの形式 | マニュアルに概要のみ | [PAO 基底](../modules/pao-basis.md)、[非局所擬ポテンシャル](../modules/nonlocal-pseudopotential.md) に記載済み |
| 二中心積分の波数格子・積分点数 | 技術ノート 15 に例示のみ | [二中心積分](../modules/two-center-overlap-kinetic.md) に記載済み |
| 軌道・原子・周期像・スピンの添字規約 | 記載なし | [疎ブロック形式](../design/sparse-block-format.md) で〔設計〕、既存実装の規約も記載 |
| 実空間グリッドの決め方（カットオフから格子数） | 技術ノート 1 に概要（係数に誤りあり） | [グリッド](../modules/grid-poisson-xc.md) に記載済み |
| k 点の生成（TRIM 点を含む） | マニュアル 15 節のみ | [k 点と Bloch 和](../modules/kpoints-bloch.md) に記載済み |

## 数値手法

| 項目 | 状況 |
|---|---|
| S の正則化（固有値の置換・次元削減）の閾値と方式 | [一般化固有値](../modules/generalized-eigen.md) に経路ごとの対応表を記載済み |
| auto-Kerker、SCF の自動調整（マニュアル 16.2、16.3） | [密度行列と電荷ミキシング](../modules/density-mixing.md) に記載済み |
| 原子密度による XC エネルギーの補正（既定で有効） | 公開資料に記載なし。[グリッド](../modules/grid-poisson-xc.md) に記載済み |
| 構造最適化の EF・RF・DIIS（マニュアル 18） | [構造最適化・MD](../modules/geometry-md.md) に記載済み（既定値はマニュアルと異なる） |
| MD の熱浴・圧力浴（NVT_NH、NPT など） | [構造最適化・MD](../modules/geometry-md.md) に記載済み |
| EXX を SCF にどう組み込むか | 既存実装では到達できない。[EXX・分散力補正](../modules/exx-vdw.md) に記載 |
| DFT-D2・D3 | [EXX・分散力補正](../modules/exx-vdw.md) に記載済み |
| O(N) 法の局所問題・化学ポテンシャル | [O(N) 法](../modules/on-methods.md) に記載済み |

## 物理機能

| 項目 | 状況 |
|---|---|
| 完全相対論・j 依存の擬ポテンシャルの扱い | [非局所擬ポテンシャル](../modules/nonlocal-pseudopotential.md) に記載済み |
| 非共線・スピン軌道相互作用でのエネルギー・力・応力 | [非共線](../modules/noncollinear-energy-force.md) に記載済み（応力は既存実装でも未対応） |
| 制約 DFT、MAE、Zeeman 項、軌道磁気モーメント | [DFT+U・磁気](../modules/dftu-magnetism.md) に記載済み |
| NEGF（非共線を含む） | [NEGF](../modules/negf.md) に記載済み |
| 電場・電荷ドープ・仮想原子・内殻準位 | [ESM・電場・内殻準位](../modules/esm-field-corehole.md) に記載済み |
| GPU 版で追加された経路 | 結果を変えないものは対象外（[並列化](../design/parallelization.md)） |

## 出力・インターフェース

| 項目 | 状況 |
|---|---|
| scfout ファイルの形式（マニュアル 70.1） | [解析機能](../modules/analysis.md) に記載済み |
| OpenMX API（マニュアル 71） | [解析機能](../modules/analysis.md) に記載済み |
| Mulliken・Voronoi・ESP 電荷、STM、仕事関数など | [解析機能](../modules/analysis.md) に記載済み |

## 未着手

| 項目 | 状況 |
|---|---|
| BoltzTraP・FermiSurfer との連携、COHP、自動 k 経路 | 未作成 |
| `Cluster-LNO`、EGAC などの一部の解法 | 一部のみ |
| 本家 OpenMX と GPU 版の差 | 未確認 |
