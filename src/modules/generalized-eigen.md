# 一般化固有値問題と S の正則化

> 状態：骨組みのみ

## 問題〔公開〕

各 k 点で一般化固有値問題 \\( H(\mathbf k)\\, c = \varepsilon\\, S(\mathbf k)\\, c \\) を解きます。
\\( H \\) は Hermite、\\( S \\) は Hermite 正定値です（数値的にはほぼ特異になりえます）。

## S の直交化と正則化〔規約〕

既存実装の主要な経路は、Cholesky 分解ではなく \\( S \\) の固有分解で直交化します。

$$
S = U \Lambda U^\dagger,\qquad X = U \tilde\Lambda^{-1/2},\qquad H' = X^\dagger H X,\qquad C = X D
$$

\\( \Lambda \\) から \\( \tilde\Lambda \\) を作る規則は経路によって異なります。

| 方式 | 規則 | 次元 |
|---|---|---|
| A | 負の固有値だけを \\( 10^{-10} \\) に置き換える | 不変 |
| B | \\( 10^{-10} \\) 未満の固有値をすべて \\( 10^{-10} \\) に置き換える | 不変 |
| C | \\( 10^{-9} \\) 未満の固有値の方向を除く（分割統治法の局所問題） | 減る |

〔未確定〕どの計算条件（Cluster/Band、共線/非共線、CPU/GPU）でどの方式が使われるかの対応表。

## 実装への要求〔設計〕

- 正則化の方式を明示的なオプションにし、既定値を文書化します。
- 置き換えが起きた場合は、残差と直交性を \\( S\_{\mathrm{eff}} = U\tilde\Lambda U^\dagger \\) に対して評価し、
  元の \\( S \\) との差も記録します。
- 必要な固有対の数 \\( n\_{\mathrm{ev}} \\) を指定できるようにします。
- SCF の反復で \\( S \\) が変わらない場合、\\( S \\) の分解を再利用できるようにします。

## テスト

- 各固有対の相対残差
  \\( \eta\_i = \lVert H c\_i - \varepsilon\_i S c\_i\rVert / \left( (\lVert H\rVert + \lvert\varepsilon\_i\rvert\\,\lVert S\rVert)\lVert c\_i\rVert \right) \\)。
- \\( \lVert C^\dagger S C - I \rVert \\)。
- 縮退がある場合は部分空間の比較。
- 良条件の \\( S \\) では LAPACK の `zhegv` 系と比較。
