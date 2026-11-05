# 実空間グリッド・Poisson・交換相関

> 状態：骨組みのみ

## 方法〔公開〕

技術ノート 1 に基づきます。

- 単位胞の格子ベクトル \\( \mathbf a\_i \\) を \\( N\_i \\) 分割した実空間格子を使います。
  \\( N\_i \\) は、各方向のカットオフエネルギー \\( E\_{\mathrm{cut}}(i) = \tfrac12 \lVert \mathbf b\_i / N\_i \rVert^2 \\) 程度
  （\\( \mathbf b\_i \\) は逆格子ベクトル）がそろい、指定値（`scf.energycutoff`、Ry）を満たすように選びます。
  マニュアル 14 節によると、通常は 150〜200 Ry で十分です。
- Hartree ポテンシャルは、原子密度の和からの差 \\( \delta n \\) について FFT で Poisson 方程式を解きます：
  \\( \delta\tilde V\_H(\mathbf q) = 4\pi\\, \delta\tilde n(\mathbf q) / \lVert\mathbf q\rVert^2 \\)。
- 交換相関ポテンシャルは格子点上で評価し、行列要素 \\( \langle\phi\_a|V|\phi\_b\rangle \\) を格子上の和で計算します。

## 〔未確定〕

- \\( N\_i \\) の決め方の正確な規則（FFT に適した数への丸めを含む）。
- \\( \mathbf q = 0 \\) 成分の扱い。
- GGA の勾配の計算方法（FFT か差分か）。
- 格子の原点と原子位置の関係。
