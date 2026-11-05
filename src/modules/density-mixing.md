# 密度行列と電荷ミキシング

> 状態：骨組みのみ

## 方法〔公開〕

- 密度行列は、固有ベクトルと占有数から \\( \rho\_{ab}(\mathbf L) = \frac{1}{N\_k}\sum\_{\mathbf k} \sum\_n f\_n(\mathbf k)\\, e^{-i\mathbf k\cdot\mathbf L}\\, c\_{an}(\mathbf k)\\, c^{\ast}\_{bn}(\mathbf k) \\) のように作ります
  （技術ノート 1 Eq. 4。位相の符号は [疎ブロック形式](../design/sparse-block-format.md) の規約に合わせて確定させる）。
- 電荷ミキシングの方式（単純、GR-Pulay、Kerker、RMM-DIIS など）は技術ノート 12 にあります。
- 電荷密度の外挿（構造最適化・MD での初期値）は技術ノート 8 にあります。

## 〔未確定〕

- 占有数のスメアリング関数（Fermi–Dirac）と化学ポテンシャルの決め方。
- auto-Kerker と SCF の自動調整の規則（マニュアル 16.2、16.3。技術ノートなし）。
- ミキシングを密度行列・実空間密度・逆空間密度のどれに適用するかの方式ごとの違い。
