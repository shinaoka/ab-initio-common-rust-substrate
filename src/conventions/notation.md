# 記法・単位・用語

## 単位〔規約〕

- 長さは bohr、エネルギーは Hartree、波数は bohr⁻¹ を内部単位とします。
- 入力ファイルでは Å・eV・Ry も使われますが、読み込み時に内部単位へ変換します。
- カットオフエネルギーを Ry で与える場合、対応する波数は \\( k = \sqrt{E\_{\mathrm{Ry}}} \\) bohr⁻¹ です
  （\\( E\_{\mathrm{Ha}} = k^2/2 \\)、\\( E\_{\mathrm{Ry}} = k^2 \\)）。

## 記号〔設計〕

| 記号 | 意味 |
|---|---|
| \\( I, J \\) | 原子（単位胞内の原子番号） |
| \\( \mathbf{R}\_I \\) | 原子 \\( I \\) の位置 |
| \\( \mathbf{L} \\) | 格子ベクトル（周期像） |
| \\( \alpha = (l, \zeta, m) \\) | 原子上の軌道。\\( l \\) は角運動量、\\( \zeta \\) は多重度（0 始まり）、\\( m \\) は磁気量子数 |
| \\( R\_{l\zeta}(r) \\) | 動径関数 |
| \\( Y\_{lm} \\) | 複素球面調和関数（Condon–Shortley 位相） |
| \\( S\_{lm} \\) | 実球面調和関数（並びは [PAO 基底](../methods/lcao-openmx/modules/pao-basis.md)） |
| \\( r\_c \\) | PAO のカットオフ半径 |
| \\( \boldsymbol\tau \\) | 二中心積分における原子間ベクトル |

## Fourier 変換の規約〔設計〕

このノートでは

$$
\tilde f(\mathbf k) = \int d^3 r\\, e^{-i\mathbf k\cdot\mathbf r} f(\mathbf r),\qquad
f(\mathbf r) = \frac{1}{(2\pi)^3}\int d^3 k\\, e^{i\mathbf k\cdot\mathbf r} \tilde f(\mathbf k)
$$

を使います。技術ノート 1 は前向き変換に \\( (2\pi)^{-3} \\) を付ける規約なので、
係数を比較する際は注意してください。動径方向の球 Bessel 変換は

$$
\tilde R\_l(k) = \int\_0^\infty dr\\, r^2 j\_l(kr) R\_l(r),\qquad
R\_l(r) = \frac{2}{\pi}\int\_0^\infty dk\\, k^2 j\_l(kr) \tilde R\_l(k)
$$

とします。
