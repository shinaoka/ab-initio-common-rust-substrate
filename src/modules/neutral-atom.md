# 中性原子ポテンシャルと射影子展開

> 状態：試作（レビュー待ち）。式は技術ノート 15（Ozaki–Kino, PRB **72**, 045121 (2005) に対応、以下 TN15）と
> 技術ノート 1（Ver. 1.3、以下 TN1）の方法を、このノートの記法で導出し直したものです。
> 数値で確かめた部分は `tools/check_neutral_atom.py` にあります。

## 記号〔設計〕

[記法](../design/notation.md) と [全エネルギーと力](energy-force.md) の記号に加えて、次を使います。

| 記号 | 意味 |
|---|---|
| \\( Z\_I \\) | 原子 \\( I \\) の種の価電子数（擬ポテンシャルの芯電荷。[非局所擬ポテンシャル](nonlocal-pseudopotential.md) 2.2 節の \\( Z\_v \\)） |
| \\( V^{\mathrm{loc}}\_I(r) \\) | 局所擬ポテンシャル。遠方で \\( -Z\_I/r \\) |
| \\( n^{a}\_I(r) \\) | 閉じ込めた孤立原子の価電子密度。\\( 4\pi\int\_0^\infty r^2 n^{a}\_I\\, dr = Z\_I \\) |
| \\( V^{a}\_{H,I}(r) \\) | \\( n^{a}\_I \\) の Hartree ポテンシャル |
| \\( V^{\mathrm{na}}\_I(r) \\) | 中性原子ポテンシャル |
| \\( r\_{c,I} \\) | PAO のカットオフ半径（PAO ファイルの `radial.cutoff.pao`） |
| \\( \boldsymbol\tau\_{IK} = \mathbf R\_K - \mathbf R\_I \\) | 原子間ベクトル（周期像の格子ベクトルを含めて考える） |
| \\( \bar R\_{K,l\zeta}(r) \\) | 原子 \\( K \\) の射影子展開の動径関数（3 節） |
| \\( c\_{K,l\zeta} \\) | その \\( V^{\mathrm{na}}\_K \\) 重みのノルム（3 節） |
| \\( \beta^{K}\_{l\zeta m}(\mathbf r) = V^{\mathrm{na}}\_K(r)\\, \bar R\_{K,l\zeta}(r)\\, S\_{lm}(\hat{\mathbf r}) \\) | 射影子（原子 \\( K \\) を原点とする） |
| \\( L\_{\max} \\)、\\( N\_{\mathrm{rad}} \\) | 射影子展開の最大角運動量と、各 \\( l \\) の動径関数の数 |

## 1. 中性原子ポテンシャル

### 1.1 定義と短距離性〔公開〕

TN1 Eq. 18、TN15（Eq. 9 の後）に従い、原子 \\( I \\) の中性原子ポテンシャルを

$$
V^{\mathrm{na}}\_I(r) = V^{\mathrm{loc}}\_I(r) + V^{a}\_{H,I}(r),\qquad
V^{a}\_{H,I}(r) = \int d^3 r'\\, \frac{n^{a}\_I(r')}{\lvert\mathbf r - \mathbf r'\rvert}
$$

と定義します。\\( n^{a}\_I \\) は PAO の生成に使ったのと同じ閉じ込めポテンシャルの下で求めた原子の価電子密度で、
\\( r > r\_{c,I} \\) で 0 です。球対称な電荷分布の Hartree ポテンシャルは

$$
V^{a}\_{H,I}(r) = \frac{4\pi}{r}\int\_0^r dr'\\, r'^2 n^{a}\_I(r') + 4\pi\int\_r^\infty dr'\\, r' n^{a}\_I(r')
$$

なので、\\( r \ge r\_{c,I} \\) では \\( V^{a}\_{H,I}(r) = Z\_I/r \\) です（Gauss の法則）。
一方、局所擬ポテンシャルは擬ポテンシャルの生成半径より外で \\( -Z\_I/r \\) に一致します。
したがって \\( V^{\mathrm{na}}\_I \\) は \\( r \ge r\_{c,I} \\) で 0 になり、カットオフ半径 \\( r\_{c,I} \\) の短距離ポテンシャルです。
物理的には、中性の原子（イオン芯と、閉じ込めた価電子雲）が外部に作るポテンシャルが 0 であることに当たります。

この性質により、\\( \sum\_I V^{\mathrm{na}}\_I \\) の行列要素は有限個の原子の組だけで決まり、
長距離の静電相互作用は差の密度 \\( \delta n = n - \sum\_I n^{a}\_I \\) の Hartree 項だけに残ります（6 節）。

### 1.2 入力データ〔規約〕

- \\( n^{a}\_I \\)：PAO ファイルの `<valence.charge.density` ブロック。各行は \\( \ln r \\)、\\( r \\)、\\( n^{a}(r) \\)（bohr⁻³）です
  （格子は動径関数の表と同じ。[PAO 基底](pao-basis.md)）。
- 読み込んだ密度は、\\( 4\pi\int r^2 n^{a}\\, dr = Z\_I \\) となるように定数倍して使います。
  配布ファイルには、表の積分が \\( Z\_I \\) と一致しないものがあります（例：C6.0 では約 4.5、\\( Z = 4 \\)）。
  このため再規格化は省略できません。
- \\( V^{\mathrm{loc}}\_I \\)：VPS ファイルの `<Pseudo.Potentials` ブロックの第 3 列（[非局所擬ポテンシャル](nonlocal-pseudopotential.md) 2.4 節）。
- \\( r\_{c,I} \\)：PAO ファイルの `radial.cutoff.pao`。基底の名前（例 `C6.0`）の数字ではありません。
- 空原子（`AtomSpecies` が 0 の種）では \\( V^{\mathrm{na}} = 0 \\)、\\( n^{a} = 0 \\) とします。

### 1.3 動径の表の作り方〔規約〕

\\( V^{\mathrm{na}}\_I \\) は VPS ファイルの動径格子 \\( \\{r\_i\\} \\) 上の表として持ちます。手順は次のとおりです。

1. 再規格化した \\( n^{a}\_I \\) から、1.1 節の式で \\( V^{a}\_{H,I}(r\_i) \\) を数値積分します
   （密度の表の外側では \\( Z\_I/r \\)）。
2. **漸近形の補正**：格子点のうち \\( r\_{c,I} \\) に最も近い点 \\( r\_{\ast} \\) で

$$
s\_I = -\frac{V^{\mathrm{loc}}\_I(r\_{\ast})}{V^{a}\_{H,I}(r\_{\ast})}
$$

   を求め、\\( V^{a}\_{H,I} \\) 全体を \\( s\_I \\) 倍します（\\( Z\_I \\) が 0 のときは補正しません）。
   1.1 節の議論から \\( s\_I \\) は数値誤差の範囲で 1 です（配布ファイルの C・O・H・Fe で \\( \lvert s\_I - 1\rvert \lesssim 10^{-4} \\)）。
   この補正で、\\( r\_{c,I} \\) で \\( V^{\mathrm{loc}}\_I + V^{a}\_{H,I} \\) が数値的にも 0 になります。
3. **滑らかな打ち切り**：

$$
V^{\mathrm{na}}\_I(r) =
\begin{cases}
\dfrac{V^{\mathrm{loc}}\_I(r) + s\_I V^{a}\_{H,I}(r)}{1 + e^{20\\,(r - r\_{c,I})}} & r \le r\_{c,I} \\\\
0 & r > r\_{c,I}
\end{cases}
$$

   （\\( r \\) は bohr）。Fermi 型の因子は \\( r\_{c,I} \\) で 1/2 ですが、2 の補正で括弧の中がほぼ 0 なので、
   \\( r\_{c,I} \\) での不連続は無視できる大きさです（配布ファイルで \\( r\_{c,I} \\) の内側 1 bohr の \\( \lvert V^{\mathrm{loc}} + sV^{a}\_H\rvert \\) は \\( 10^{-4} \\) Hartree 程度以下）。
   \\( r > r\_{c,I} \\) での値は補間に関係なく厳密に 0 として扱います。

\\( V^{a}\_{H,I} \\) 自体（補正後）も、[全エネルギーと力](energy-force.md) の遮蔽された原子間項で使います。

〔設計〕補間の方式は定めません。\\( r < r\_1 \\)（最初の格子点より内側）では \\( V^{\mathrm{na}} \\) を \\( a r^2 + b \\) の形で
滑らかに延長します（原点で正則なので）。既存実装もこの形の延長を使っています〔規約〕。

### 1.4 波数空間の表〔規約〕

4 節と 3.5 節の波数空間の積分には、\\( V^{\mathrm{na}}\_I \\) と \\( \beta^{K} \\) の動径関数の球 Bessel 変換

$$
\tilde V^{\mathrm{na}}\_I(k) = \int\_0^{\infty} dr\\, r^2 j\_0(kr)\\, V^{\mathrm{na}}\_I(r)
$$

を使います（[記法](../design/notation.md) の規約）。既存実装は、実空間の積分を \\( r\_{c,I} + 0.5 \\) bohr で打ち切り、
波数は \\( 10^{-6} \le k \le k\_{\max} \\)（\\( k\_{\max} = \sqrt{E} \\)、\\( E \\) は `1DFFT.EnergyCutoff`、既定 3600 Ry で 60 bohr⁻¹）の
区間の Gauss–Legendre 300 点の上で表にします。波数の積分も同じ 300 点で行います。
[二中心積分](two-center-overlap-kinetic.md) の重なり・運動エネルギーの等間隔格子とは別の格子です。

## 2. 行列要素の分類〔公開〕〔規約〕

ハミルトニアンの中性原子項は

$$
V^{\mathrm{na}}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{K} \langle \phi\_{I\alpha} \vert V^{\mathrm{na}}\_K \vert \phi\_{J\beta}(\cdot - \boldsymbol\tau\_{IJ\mathbf L}) \rangle
$$

で、和は \\( V^{\mathrm{na}}\_K \\) の台が両方の基底関数と重なるすべての原子 \\( K \\)（周期像を含む）にわたります。
\\( K \\) と \\( I \\)、\\( J \\) の関係で、次のように計算方法を分けます。

| 場合 | 例 | 計算方法 | 出典 |
|---|---|---|---|
| 同じ原子上の基底対（\\( J = I \\)、\\( \mathbf L = 0 \\)） | \\( K = I \\)（一中心）、\\( K \ne I \\)（二中心） | 波数空間で直接（4 節） | TN15、〔規約〕 |
| 異なる原子の基底対で、\\( K = I \\) または \\( K = J \\) | 二中心 | 射影子展開（3 節）。基底関数が展開の動径関数の張る空間に入るので厳密 | TN15 Eq. 16、〔規約〕 |
| 異なる原子の基底対で、\\( K \ne I, J \\) | 三中心 | 射影子展開（近似） | TN15 Eq. 13 |
| `scf.ProExpn.VNA off` のとき | すべて | 実空間格子上の和（5 節） | マニュアル、〔規約〕 |

〔公開〕TN15 は、\\( K = I \\) または \\( K = J \\) の二中心の場合に、展開が \\( V^{\mathrm{na}}\_K\vert\phi\rangle \\) を厳密に再現する性質（TN15 Eq. 16）を使うこと、
\\( I = J \ne K \\) の場合は \\( V^{\mathrm{na}}\_K \\) と動径関数の積を波数空間で扱って別に計算することを述べています。

〔規約〕既存実装（`scf.ProExpn.VNA on`、既定）の具体的な分け方は次のとおりです。

- **同じ原子上のブロック**（\\( J = I \\)、\\( \mathbf L = 0 \\)）：\\( K = I \\) を含むすべての \\( K \\) について、4 節の直接計算の和を使います。
  射影子展開の値は使いません。
- **それ以外のブロック**（\\( J \ne I \\) または \\( \mathbf L \ne 0 \\)。同じ原子の周期像どうしも含む）：
  \\( K = I \\)、\\( K = J \\) を含むすべての \\( K \\) について射影子展開の和をとり、3.6 節の滑らかな打ち切り因子を掛けます。
  \\( K = I, J \\) の項を別の式で計算することはしません（展開が厳密なので、結果は同じです）。

## 3. 射影子展開

### 3.1 展開の形〔公開〕

TN15 Eq. 13–15（TN1 Eq. 32–34）に従い、原子 \\( K \\) の \\( V^{\mathrm{na}}\_K \\) を分離型の演算子

$$
\hat V^{\mathrm{na}}\_K = \sum\_{l=0}^{L\_{\max}} \sum\_{m=-l}^{l} \sum\_{\zeta=1}^{N\_{\mathrm{rad}}}
\frac{\lvert \beta^{K}\_{l\zeta m} \rangle \langle \beta^{K}\_{l\zeta m} \rvert}{c\_{K,l\zeta}},\qquad
\beta^{K}\_{l\zeta m}(\mathbf r) = V^{\mathrm{na}}\_K(r)\\, \bar R\_{K,l\zeta}(r)\\, S\_{lm}(\hat{\mathbf r})
$$

で近似します。動径関数 \\( \bar R\_{l\zeta} \\)（以下 \\( K \\) を省略）は、種の関数 \\( \chi\_{l\zeta} \\)（3.2 節）を
\\( V^{\mathrm{na}}\_K \\) を重みとする内積

$$
(f, g)\_V = \int\_0^\infty dr\\, r^2 f(r)\\, V^{\mathrm{na}}\_K(r)\\, g(r)
$$

で Gram–Schmidt 直交化したものです：

$$
\bar R\_{l\zeta} = \chi\_{l\zeta} - \sum\_{\eta < \zeta} \frac{(\bar R\_{l\eta}, \chi\_{l\zeta})\_V}{c\_{l\eta}}\\, \bar R\_{l\eta},\qquad
c\_{l\zeta} = (\bar R\_{l\zeta}, \bar R\_{l\zeta})\_V .
$$

**性質（TN15 Eq. 16）**：動径関数が \\( \mathrm{span}\\{\bar R\_{l\zeta}\\}\_\zeta \\) に入る関数
\\( f(\mathbf r) = R(r) S\_{lm}(\hat{\mathbf r}) \\)（\\( l \le L\_{\max} \\)）については \\( \hat V^{\mathrm{na}}\_K f = V^{\mathrm{na}}\_K f \\) が厳密に成り立ちます。
実際、\\( R = \sum\_\zeta a\_\zeta \bar R\_{l\zeta} \\) と書くと、直交性 \\( (\bar R\_{l\zeta'}, \bar R\_{l\zeta})\_V = \delta\_{\zeta\zeta'} c\_{l\zeta} \\) と
球面調和関数の直交性から \\( \langle \beta^{K}\_{l'\zeta' m'} \vert f \rangle = \delta\_{ll'}\delta\_{mm'} a\_{\zeta'} c\_{l\zeta'} \\) となり、

$$
\hat V^{\mathrm{na}}\_K f = \sum\_{\zeta'} V^{\mathrm{na}}\_K \bar R\_{l\zeta'} S\_{lm}\\, a\_{\zeta'} = V^{\mathrm{na}}\_K f
$$

です。したがって、原子 \\( K \\) 上の基底関数が展開の空間に入っていれば、任意の \\( \phi \\) について
\\( \langle \phi \vert \hat V^{\mathrm{na}}\_K \vert \phi\_{K\alpha} \rangle = \langle \phi \vert V^{\mathrm{na}}\_K \vert \phi\_{K\alpha} \rangle \\) です。
三中心の項は \\( L\_{\max} \\)、\\( N\_{\mathrm{rad}} \\) を増やすと厳密な値に収束します（TN15 Fig. 2）。

- \\( V^{\mathrm{na}} \\) は符号が一定とは限らない（多くの原子でほぼ負）ので、\\( (\cdot,\cdot)\_V \\) は内積としては不定値です。
  \\( c\_{l\zeta} \\) は負にも 0 に近くにもなりえます。この展開は \\( V^{\mathrm{na}}\_K \\) の符号を保つ必要がない点に注意してください。
- 分離型なので、三中心の行列要素は二中心積分 \\( \langle \phi \vert \beta \rangle \\) の積の和になり、
  [非局所擬ポテンシャル](nonlocal-pseudopotential.md) と同じ組み立て方が使えます（射影子エネルギー \\( D \\) が対角の \\( 1/c \\) に当たる）。
- \\( \hat V^{\mathrm{na}}\_K \\) は回転不変です（各 \\( l \\) で \\( m \\) について和をとるので）。\\( S\_{lm} \\) を複素の \\( Y\_{lm} \\) に替えても同じ演算子です。

### 3.2 種の関数〔公開〕〔規約〕

〔公開〕TN15 は、種の関数として、\\( n^{a} \\) と同じ閉じ込めポテンシャルの下で求めた基底状態と励起状態の擬原子軌道を使うと述べています。

〔規約〕既存実装の種の関数は次のとおりです。PAO ファイルに格納された動径関数を \\( R^{\mathrm{PAO}}\_{l\zeta} \\)
（\\( l = 0,\dots,L\_p \\)、\\( \zeta = 0,\dots,M\_p - 1 \\)、\\( L\_p \\) は `PAO.Lmax`、\\( M\_p \\) は `PAO.Mul`）、
\\( u(r) = V^{\mathrm{na}}\_K(r)/10 + 10^{-13} \\) とし、\\( \zeta \\) を 0 始まりで数えます。

| 条件 | 種の関数 \\( \chi\_{l\zeta}(r) \\) |
|---|---|
| \\( l \le L\_p \\)、\\( \zeta < M\_p \\) | \\( R^{\mathrm{PAO}}\_{l\zeta}(r) \\) |
| \\( l \le L\_p \\)、\\( \zeta \ge M\_p \\) | \\( R^{\mathrm{PAO}}\_{l0}(r)\\, u(r)^{\zeta} \\) |
| \\( l > L\_p \\)、\\( \zeta < M\_p \\) | \\( R^{\mathrm{PAO}}\_{L\_p\zeta}(r)\\, r^{\\,l - L\_p} \\) |
| \\( l > L\_p \\)、\\( \zeta \ge M\_p \\) | \\( R^{\mathrm{PAO}}\_{L\_p, M\_p-1}(r)\\, u(r)^{\\,\zeta - M\_p + 1} \\) |

- 使う PAO は、計算で選んだ基底（例 `s2p2d1`）ではなく、ファイルに格納されたものすべての先頭からです。
  配布されている PAO ファイルは \\( M\_p = 15 \\)（1 つだけ 20、1 つだけ 5）なので、既定の \\( N\_{\mathrm{rad}} = 12 \\) では
  \\( l \le L\_p \\) の種はすべて PAO そのものです。
- 基底関数は PAO ファイルの先頭の動径関数（またはその線形結合）なので、各 \\( l \le \min(L\_p, L\_{\max}) \\) で選んだ動径関数の数が
  \\( N\_{\mathrm{rad}} \\) 以下なら、3.1 節の性質により二中心の項は厳密です。
- \\( l > L\_p \\) の種は、\\( l = L\_p \\) の PAO に \\( r^{l-L\_p} \\) を掛けて原点付近の振る舞い \\( r^l \\) を合わせたものです。
  既定の \\( L\_{\max} \\)（3.4 節）は通常 \\( L\_p \\) より大きいので、この行は既定の計算で実際に使われます。
- 動径関数は VPS の動径格子に補間して扱います（[PAO 基底](pao-basis.md) の範囲外の値の規則に従う）。

### 3.3 直交化の数値的な扱い〔規約〕〔設計〕

〔規約〕既存実装の手順の要点は次のとおりです。

- 種の関数は、直交化の前に通常のノルム \\( \int r^2 \chi^2 dr \\) で規格化します。
- 直交化は 3.1 節の式そのもの（係数に元の種 \\( \chi\_{l\zeta} \\) を使う古典的 Gram–Schmidt）です。
- \\( \lvert c\_{l\zeta}\rvert < 10^{-15} \\) の項は \\( 1/c\_{l\zeta} = 0 \\) として捨てます。空原子ではすべての項を捨てます。
- 動径積分は、VPS の動径格子の最初の点から最後の点までの区間を Gauss–Legendre 300 点で行います。
- 規格化の取り方（\\( \bar R \\) の定数倍）は \\( \lvert\beta\rangle\langle\beta\rvert/c \\) で打ち消すので、展開に影響しません。

〔設計〕再実装では、直交化を修正 Gram–Schmidt や、\\( (\cdot,\cdot)\_V \\) の Gram 行列の対角化で行ってもよいです
（厳密な算術では同じ展開になります）。ただし、種の関数がほぼ一次従属になる（\\( c \\) が非常に小さい）と
\\( 1/c \\) が大きくなって丸め誤差が増えるので、上の閾値のような安全策を必ず入れます。

### 3.4 打ち切りの値〔公開〕〔規約〕

〔公開〕TN15 は、H₂O 分子での収束の結果から、\\( L\_{\max} \\) を「基底の最大の角運動量 + 4」、\\( N\_{\mathrm{rad}} = 4 \\) とすることを推奨しています。
水分子では \\( L\_{\max} \approx 6 \\)、\\( N\_{\mathrm{rad}} \approx 4 \\) で収束し、どちらも最低 2 は必要だとしています。

〔規約〕既存実装の既定値はこれより大きく、入力キーワードで変えられます。

| キーワード | 既定値 | 意味 |
|---|---|---|
| `scf.ProExpn.VNA` | on | 射影子展開を使うか。off なら実空間格子で計算（マニュアル 8.2 節に記載） |
| `scf.BufferL.VNA` | 6 | \\( L\_{\max} = \max\_{\text{全種}}(\text{基底の最大の } l) + \\) この値 |
| `scf.RadialF.VNA` | 12 | \\( N\_{\mathrm{rad}} \\) |

- \\( L\_{\max} \\)、\\( N\_{\mathrm{rad}} \\) は全種で共通です（種ごとではなく、計算に現れるすべての種の基底の最大 \\( l \\) で決まる）。
  例：基底が `s2p2d1` までなら \\( L\_{\max} = 8 \\)。
- 後の 2 つのキーワードはマニュアルのキーワード一覧には見当たりません。

〔設計〕再実装は既存実装の既定値（6 と 12）に合わせます。fixture との照合も既定値で行います。

### 3.5 射影子との重なり〔公開〕

\\( \langle \phi\_{I\alpha} \vert \beta^{K}\_{l\zeta m} \rangle \\) は、[二中心積分](two-center-overlap-kinetic.md) の重なりの式で、
\\( \phi\_b \\) の動径関数を \\( V^{\mathrm{na}}\_K \bar R\_{K,l\zeta} \\) に、原子間ベクトルを \\( \boldsymbol\tau\_{IK} \\) に替えたものです（TN1 の方法）。
\\( V^{\mathrm{na}}\_K \bar R\_{K,l\zeta} \\) の球 Bessel 変換を 1.4 節と同じ格子で表にしておきます。
\\( \lVert\boldsymbol\tau\_{IK}\rVert \ge r\_{c,I} + r\_{c,K} \\) では厳密に 0 です。

### 3.6 三中心の和と滑らかな打ち切り〔規約〕

〔規約〕異なる原子のブロックは

$$
V^{\mathrm{na}}\_{I\alpha,J\beta}(\mathbf L) = d\big(\lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert;\\, r\_{c,I} + r\_{c,J}\big)
\sum\_{K} \sum\_{l\zeta m} \frac{\langle \phi\_{I\alpha} \vert \beta^{K}\_{l\zeta m} \rangle \langle \beta^{K}\_{l\zeta m} \vert \phi\_{J\beta}(\cdot-\boldsymbol\tau\_{IJ\mathbf L}) \rangle}{c\_{K,l\zeta}}
$$

です。\\( K \\)（周期像を含む）は、\\( I \\) とも \\( J \\) とも原子対の範囲（[疎ブロック形式](../design/sparse-block-format.md)、
距離が \\( r\_c \\) の和以下）にあるものにわたります。この範囲の外では片方の重なりが 0 なので、和の打ち切りは厳密です。

打ち切り因子 \\( d \\) は、バッファ幅 1 bohr の \\( C^2 \\) の階段関数です。\\( x = R - R\_c \\) として

$$
d(R; R\_c) =
\begin{cases}
1 & R \le R\_c - 1 \\\\
{}-x^3\\,(10 + 15x + 6x^2) & R\_c - 1 < R \le R\_c \\\\
0 & R > R\_c
\end{cases}
$$

（\\( R \\) は bohr）。\\( x = -1 \\) で値 1・1 階と 2 階の微分 0、\\( x = 0 \\) で値・1 階・2 階の微分 0 になる 5 次多項式です。

- **理由**：厳密な行列要素 \\( \langle\phi\_I\vert V^{\mathrm{na}}\_K\vert\phi\_J\rangle \\) は、\\( \phi\_I \\) と \\( \phi\_J \\) の台が重ならない
  \\( \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \ge r\_{c,I} + r\_{c,J} \\) で 0 ですが、分離型の近似ではこの性質が失われ、
  原子対の範囲の端で値が不連続に 0 になります。因子 \\( d \\) は、対が範囲に出入りするときにエネルギーと力を連続にします。
- 同じ原子上のブロック（2 節の直接計算）には掛けません（\\( d(0) = 1 \\) でもあります）。
- 因子は行列要素そのものに掛かるので、ハミルトニアン、`Una`、力のすべてに一貫して現れます。力ではその微分も必要です（7 節）。

## 4. 同じ原子上の基底対と別原子の \\( V^{\mathrm{na}} \\)〔公開〕〔規約〕

TN15 が波数空間で扱うとしている場合（\\( I = J \\)、\\( \mathbf L = 0 \\)、\\( K \\) は任意）の式を、このノートの規約で導きます。

原子 \\( I \\) を原点とし、\\( \boldsymbol\tau = \boldsymbol\tau\_{IK} \\)、
\\( \phi\_a = R\_a(r) Y\_{lm} \\)、\\( \phi\_b = R\_b(r) Y\_{l'm'} \\) とします。角度部分の積を

$$
Y^{\ast}\_{lm} Y\_{l'm'} = \sum\_{LM} G(lm;\\, l'm',\\, LM)\\, Y^{\ast}\_{LM}
$$

と展開し（\\( G \\) は [二中心積分](two-center-overlap-kinetic.md) の Gaunt 係数、和は \\( M = m - m' \\) だけ）、
\\( V^{\mathrm{na}}\_K(\lvert\mathbf r - \boldsymbol\tau\rvert) \\) を Fourier 積分で書いて平面波展開を 2 回使うと、

$$
\langle \phi\_a \vert V^{\mathrm{na}}\_K(\cdot - \boldsymbol\tau) \vert \phi\_b \rangle
{}= 8 \sum\_{LM} G(lm;\\, l'm',\\, LM)\\, Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})
\int\_0^\infty dk\\, k^2 j\_L(k\tau)\\, \tilde P^{ab}\_L(k)\\, \tilde V^{\mathrm{na}}\_K(k)
$$

$$
\tilde P^{ab}\_L(k) = \int\_0^\infty dr\\, r^2 j\_L(kr)\\, R\_a(r) R\_b(r)
$$

が得られます。導出の要点：\\( \int d^3r\\, f(r) Y^{\ast}\_{LM}(\hat{\mathbf r}) e^{i\mathbf k\cdot\mathbf r} = 4\pi i^L Y^{\ast}\_{LM}(\hat{\mathbf k}) \tilde f\_L(k) \\)、
\\( e^{-i\mathbf k\cdot\boldsymbol\tau} \\) の展開の係数 \\( (-i)^L \\) と打ち消し合い位相が残りません。
係数 8 は \\( (4\pi)^3/(2\pi)^3 \\) です（3 次元の Fourier 変換 \\( 4\pi\tilde V^{\mathrm{na}}\_K \\) の \\( 4\pi \\) を含む）。

- \\( \tau = 0 \\)（\\( K = I \\)、一中心）：\\( j\_L(0) = \delta\_{L0} \\) なので
  \\( \langle\phi\_a\vert V^{\mathrm{na}}\_I\vert\phi\_b\rangle = \delta\_{ll'}\delta\_{mm'}\int\_0^\infty dr\\, r^2 R\_a R\_b V^{\mathrm{na}}\_I \\) です。
  方向 \\( \hat{\boldsymbol\tau} \\) は不定なので特別に扱います。
- 実球面調和関数の基底では、[二中心積分](two-center-overlap-kinetic.md) 4 節と同じ変換で実数にします。
- \\( \tau \ge r\_{c,I} + r\_{c,K} \\) では厳密に 0 です。
- 〔規約〕既存実装は \\( \tilde P^{ab}\_L \\) を \\( \tilde V^{\mathrm{na}}\_K \\) と同じ Gauss–Legendre の波数格子で表にし、同じ格子で \\( k \\) 積分します。

> **検算**：上の式は、Gauss 型の動径関数と解析的なモデル \\( V^{\mathrm{na}} \\)（下のテストの節）について、
> \\( (l, l') \\) の 4 組で三次元格子上の直接積分と 6 桁以上一致します（`tools/check_neutral_atom.py` の (2)）。

## 5. 格子での評価（`scf.ProExpn.VNA off`）〔公開〕〔規約〕

〔公開〕マニュアル（8.2 節）によると、`scf.ProExpn.VNA` を off にすると、\\( V^{\mathrm{na}} \\) の行列要素を実空間格子上で計算します。

〔規約〕このとき \\( \sum\_K V^{\mathrm{na}}\_K(\lvert\mathbf r\_p - \mathbf R\_K\rvert) \\) を格子点上で作り、
\\( \delta V\_H + v\_{\mathrm{xc},\sigma} \\) に加えた有効ポテンシャルとして [グリッド](grid-poisson-xc.md) の 5 節の方法で行列要素にします。
`Una` は \\( \Delta V\sum\_p n(\mathbf r\_p) V^{\mathrm{na}}(\mathbf r\_p) \\) です（[全エネルギーと力](energy-force.md)）。
既定の on のときは、\\( V^{\mathrm{na}} \\) を格子には載せません。

TN15（Fig. 3、Table I）は、格子で評価する方法では格子に対する原子の位置によってエネルギーが振動し（egg-box 効果）、
射影子展開でそれが大きく減ることを示しています。

## 6. 全エネルギーの分割での役割〔公開〕

全エネルギーの分割の導出は [全エネルギーと力](energy-force.md) にあります。ここでは中性原子ポテンシャルに関わる点だけをまとめます（TN1 Eq. 18–23）。

- \\( E\_{\mathrm{na}} = \int n \sum\_K V^{\mathrm{na}}\_K \\) は、射影子展開のときは密度行列と 2–3 節の行列要素の積和（`Una`）です。
- 静電エネルギーのうち長距離の部分は、差の密度 \\( \delta n = n - \sum\_I n^{a}\_I \\) とその Hartree ポテンシャル \\( \delta V\_H \\) の
  \\( E\_{\delta\mathrm{ee}} = \frac12\int\delta n\\,\delta V\_H \\) だけに残り、これを FFT で解きます（[グリッド](grid-poisson-xc.md)）。
  \\( n^{a}\_I \\) は 1.2 節で再規格化した同じ密度を使います。
- イオン芯どうしの反発と原子密度どうしの Hartree 項は、遮蔽された短距離の対ポテンシャル

$$
E\_{\mathrm{scc}} = \frac12\sum\_{I,J}\left[\frac{Z\_I Z\_J}{\lVert\mathbf R\_I - \mathbf R\_J\rVert}(1 - \delta\_{IJ}) - \int d^3r\\, n^{a}\_I\\, V^{a}\_{H,J}\right]
$$

  にまとまります（TN1 Eq. 22。\\( J \\) は周期像を含み、\\( I = J \\) は同じ像の場合）。
  各対の項は原子密度の台が重ならない距離で 0 です。既存実装での評価（`Ucore` と `UH0`、端点での規格化を含む）は
  [全エネルギーと力](energy-force.md) の「遮蔽された原子間項」の節にあります。

〔設計〕球対称な 2 つの電荷の Coulomb 積分は、1 次元の積分

$$
\int d^3r\\, n^{a}\_I(\mathbf r)\\, V^{a}\_{H,J}(\mathbf r - \mathbf d) = \frac{2}{\pi}\int\_0^\infty dk\\, j\_0(kd)\\, N\_I(k) N\_J(k),\qquad
N\_I(k) = 4\pi\int\_0^\infty dr\\, r^2 j\_0(kr)\\, n^{a}\_I(r)
$$

でも計算できます（\\( N\_I(0) = Z\_I \\) なので、遠方で \\( Z\_I Z\_J/d \\)）。距離の表を作る場合の代替として使えます。
Gauss 型の電荷について解析解と \\( 10^{-14} \\) で一致することを確認しています（`tools/check_neutral_atom.py` の (1)）。

## 7. 力〔公開〕〔規約〕

全体の形（密度行列を固定した行列要素の微分と、エネルギー重み付き密度行列による Pulay の項）は
[全エネルギーと力](energy-force.md) の「力」の節にあります（TN1 Eq. 59–61、75）。ここでは中性原子の行列要素の微分を与えます。
\\( g\_s \\) と密度行列の記号は同章に従います。

### 7.1 射影子展開のブロック

\\( A^{K}\_{I\alpha,p} = \langle \phi\_{I\alpha} \vert \beta^{K}\_p \rangle \\)（\\( p = (l\zeta m) \\)）は \\( \boldsymbol\tau\_{IK} \\) だけの関数なので、

$$
\frac{\partial A^{K}\_{I\alpha,p}}{\partial \mathbf R\_K} = \nabla\_{\boldsymbol\tau} A^{K}\_{I\alpha,p}(\boldsymbol\tau\_{IK}),\qquad
\frac{\partial A^{K}\_{I\alpha,p}}{\partial \mathbf R\_I} = -\nabla\_{\boldsymbol\tau} A^{K}\_{I\alpha,p}(\boldsymbol\tau\_{IK})
$$

です（\\( \nabla\_{\boldsymbol\tau} \\) は [二中心積分](two-center-overlap-kinetic.md) 6 節の解析微分。\\( K = I \\) のときは 2 つが打ち消して 0）。
3.6 節の行列要素を \\( V = d\\, M \\)、\\( M = \sum\_K\sum\_p A^{K}\_{I\alpha,p} c\_p^{-1} A^{K}\_{J\beta,p} \\) と書くと、原子 \\( X \\) についての微分は

$$
\frac{\partial V^{\mathrm{na}}\_{I\alpha,J\beta}}{\partial\mathbf R\_X}
{}= d \sum\_K\sum\_p \frac{1}{c\_p}\left[\frac{\partial A^{K}\_{I\alpha,p}}{\partial\mathbf R\_X} A^{K}\_{J\beta,p} + A^{K}\_{I\alpha,p}\frac{\partial A^{K}\_{J\beta,p}}{\partial\mathbf R\_X}\right]
{}+ M\\, d'(\tau\_{IJ})\\, \frac{\partial \tau\_{IJ}}{\partial \mathbf R\_X}
$$

です。\\( \tau\_{IJ} = \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \\)、\\( \partial\tau\_{IJ}/\partial\mathbf R\_J = \hat{\boldsymbol\tau}\_{IJ\mathbf L} = -\partial\tau\_{IJ}/\partial\mathbf R\_I \\)、
\\( d'(R) = -30x^2(1+x)^2 \\)（\\( R\_c - 1 < R \le R\_c \\)、それ以外は 0）です。
エネルギーへの寄与は、これに \\( g\_s\sum\_\sigma\rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\) を掛けた和です。
〔規約〕既存実装は打ち切り因子の微分の項を含めます。

### 7.2 同じ原子上のブロック

\\( W^{IK}\_{\alpha\beta}(\boldsymbol\tau\_{IK}) = \langle \phi\_{I\alpha} \vert V^{\mathrm{na}}\_K \vert \phi\_{I\beta} \rangle \\)（4 節の式）は \\( \boldsymbol\tau\_{IK} \\) だけの関数です。

$$
\frac{\partial W^{IK}}{\partial \mathbf R\_K} = \nabla\_{\boldsymbol\tau} W^{IK},\qquad
\frac{\partial W^{IK}}{\partial \mathbf R\_I} = -\nabla\_{\boldsymbol\tau} W^{IK}\quad(K \ne I)
$$

で、\\( K = I \\) の項は配置によらない定数です。\\( \nabla\_{\boldsymbol\tau} \\) は 4 節の式の各項
\\( f\_L(\tau) Y^{\ast}\_{LM}(\hat{\boldsymbol\tau}) \\)（\\( f\_L(\tau) = \int dk\\, k^2 j\_L(k\tau)\tilde P^{ab}\_L\tilde V^{\mathrm{na}}\_K \\)）に
[二中心積分](two-center-overlap-kinetic.md) 6 節の勾配の公式をそのまま当てはめ、
\\( f\_L'(\tau) = \int dk\\, k^3 j\_L'(k\tau)\tilde P^{ab}\_L\tilde V^{\mathrm{na}}\_K \\) を使います。
〔規約〕既存実装もこの解析微分を使い、\\( \sin\theta \\) が小さいときは \\( \varphi \\) 方向の項を落とします。

これで [全エネルギーと力](energy-force.md) の〔未確定〕にある「射影子展開を使わない組み合わせの行列要素の微分」は、
既定の `scf.ProExpn.VNA on` については上の式で定まります。

### 7.3 その他

- 格子で評価する場合（5 節）の力は [全エネルギーと力](energy-force.md) の格子上の項（第 4 項）です。
- \\( E\_{\mathrm{scc}} \\) の力は [全エネルギーと力](energy-force.md) の「遮蔽された原子間項の微分」です。
- 並進不変性から、各ブロックについて関係する原子の微分の和は 0 です（7.1 節では \\( X = I, J, K \\) の和、7.2 節では \\( X = I, K \\) の和）。

## 8. 数値パラメータのまとめ〔規約〕

| 項目 | 値 |
|---|---|
| \\( V^{\mathrm{na}} \\) の打ち切り | \\( r\_c \\)（`radial.cutoff.pao`）で 0。Fermi 型因子の急峻さ 20 bohr⁻¹ |
| 漸近形の補正 | \\( r\_c \\) に最も近い VPS 格子点で \\( V^{\mathrm{loc}} + sV^{a}\_H = 0 \\) |
| 射影子展開 | \\( L\_{\max} = \max(l\_{\mathrm{basis}}) + 6 \\)、\\( N\_{\mathrm{rad}} = 12 \\)（`scf.BufferL.VNA`、`scf.RadialF.VNA`） |
| 直交化の動径積分 | VPS 格子の全範囲、Gauss–Legendre 300 点 |
| 小さい \\( c \\) の切り捨て | \\( \lvert c\rvert < 10^{-15} \\) |
| 球 Bessel 変換 | 実空間 \\( r \le r\_c + 0.5 \\) bohr、波数 \\( [10^{-6}, \sqrt{E\_{\mathrm{1DFFT}}}] \\) を Gauss–Legendre 300 点 |
| 打ち切り因子 | \\( C^2 \\) の 5 次多項式、バッファ 1 bohr、\\( R\_c = r\_{c,I} + r\_{c,J} \\) |

## テスト

### 解析的なモデル〔設計〕

Gauss 型の原子密度と、誤差関数型の局所ポテンシャル

$$
n^{a}(r) = Z\left(\frac{\alpha}{\pi}\right)^{3/2} e^{-\alpha r^2},\qquad
V^{\mathrm{loc}}(r) = -Z\frac{\mathrm{erf}(\sqrt{\beta}\\,r)}{r}\quad(\beta > \alpha)
$$

から \\( V^{\mathrm{na}}(r) = Z\\,[\mathrm{erf}(\sqrt\alpha\\, r) - \mathrm{erf}(\sqrt\beta\\, r)]/r \\) が得られます。
有限の台は持ちませんが、Gauss 型で減衰するので、1.3 節の打ち切りを除く部分の検査に使えます。

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 1.1 節の動径積分による \\( V^{a}\_H \\) とモデルの \\( Z\\,\mathrm{erf}(\sqrt\alpha r)/r \\) | \\( 10^{-10} \\) |
| 解析 | 6 節の球対称電荷の Coulomb 積分と \\( Z\_1 Z\_2\\,\mathrm{erf}(\sqrt{\mu}d)/d \\)（\\( \mu = a\_1 a\_2/(a\_1+a\_2) \\)） | \\( 10^{-12} \\) |
| 解析 | 4 節の式と三次元格子上の直接積分（Gauss 型の \\( R\_a, R\_b \\)、\\( l, l' \le 2 \\)） | \\( 10^{-6} \\)（格子の精度で決まる） |
| 解析 | 4 節の式の \\( \tau \to 0 \\) の極限と \\( \int r^2 R\_a R\_b V^{\mathrm{na}} \\) | \\( 10^{-10} \\) |
| 解析 | 打ち切り因子：\\( x = -1, 0 \\) での値・1 階・2 階微分 | 丸め誤差 |

### 射影子展開〔設計〕

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 不変量 | \\( (\bar R\_{l\zeta}, \bar R\_{l\zeta'})\_V = \delta\_{\zeta\zeta'} c\_{l\zeta} \\) | \\( 10^{-12} \\)（相対） |
| 不変量 | 種の関数 \\( \chi \\) について \\( \hat V^{\mathrm{na}}\chi = V^{\mathrm{na}}\chi \\)：動径の内積 \\( \sum\_\zeta (\bar R\_\zeta, \chi)\_V (\bar R\_\zeta, g)\_V / c\_\zeta = (\chi, g)\_V \\) を任意の \\( g \\) で | \\( 10^{-10} \\)（相対） |
| 不変量 | 二中心（\\( K = I \\)）の展開の値と、直接の三次元積分 | 格子の精度。`tools/check_neutral_atom.py` の (3) では 6 桁 |
| 収束 | 三中心の値が \\( L\_{\max} \\)、\\( N\_{\mathrm{rad}} \\) とともに直接積分に近づく | 傾向を確認（下の例） |
| 不変量 | 行列要素の対称性 \\( V^{\mathrm{na}}\_{I\alpha,J\beta}(\mathbf L) = V^{\mathrm{na}}\_{J\beta,I\alpha}(-\mathbf L) \\) | 丸め誤差 |
| 不変量 | 構造の回転で行列が実球面調和関数の回転行列で変換される | \\( 10^{-10} \\) |
| 不変量 | 原子対の範囲の端（\\( \tau = r\_{c,I} + r\_{c,J} \\)）で異なる原子のブロックが 0、その内側 1 bohr で滑らか | 厳密 |

収束の例（`tools/check_neutral_atom.py` の (3)。s 型と p 型の Gauss 関数、\\( V^{\mathrm{na}} \\) は第 3 の中心）：
\\( N\_{\mathrm{rad}} = 4 \\) で \\( L\_{\max} = 0, 2, 4, 6 \\) のとき、直接積分との相対誤差はおよそ
\\( 6\times10^{-2} \\)、\\( 6\times10^{-3} \\)、\\( 2\times10^{-4} \\)、\\( 1\times10^{-5} \\) です。
\\( N\_{\mathrm{rad}} \\) を 2 から 6 に増やす効果は \\( 10^{-3} \\) 程度で、\\( L\_{\max} \\) の効果の方が大きいことは TN15 Fig. 2 の傾向と合います。

### 実際のデータでの検査〔設計〕

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 不変量 | 再規格化後の \\( 4\pi\int r^2 n^{a} = Z \\) | \\( 10^{-10} \\) |
| 不変量 | 1.3 節の \\( s \\) が 1 に近い（C、O、H、Fe などの配布ファイル） | \\( \lvert s-1\rvert < 10^{-3} \\) |
| 不変量 | \\( r\_c \\) の内側 1 bohr で \\( \lvert V^{\mathrm{loc}} + sV^{a}\_H\rvert \\) が小さい | \\( 10^{-3} \\) Hartree |
| 有限差分 | \\( \rho \\) を固定した `Una` の原子変位による中心差分と、7 節の解析的な力（打ち切り因子の微分を含む） | \\( 10^{-7} \\) Hartree/bohr |
| 有限差分 | 原子対の範囲の端をまたぐ変位でも上の一致が保たれる | 同上 |
| 不変量 | `Una` の力の寄与の総和が 0 | \\( 10^{-10} \\) |

### fixture との照合〔設計〕

- 既存実装の \\( V^{\mathrm{na}}\_I(r) \\)、\\( \bar R\_{l\zeta} \\)、\\( c\_{l\zeta} \\)（種ごと）を書き出して比べます。
  \\( \bar R \\) は規格化の取り方が実装で異なりうるので、\\( \beta\beta/c \\) の組み合わせ（例：動径の核 \\( \sum\_\zeta V\bar R\_\zeta(r)V\bar R\_\zeta(r')/c\_\zeta \\)）で比べます。
- 行列要素：H₂O 分子（`H6.0-s1p1`、`O6.0-s1p1` など）と Si 結晶で、同じ原子上のブロックと異なる原子のブロックの \\( V^{\mathrm{na}} \\) を比べます。
- エネルギー：`Una`（`scf.ProExpn.VNA` の on と off の両方）。
- 許容誤差：\\( V^{\mathrm{na}}(r) \\) は \\( 10^{-8} \\) Hartree、行列要素は \\( 10^{-7} \\) Hartree、`Una` は \\( 10^{-6} \\) Hartree を目安とします
  （同じ \\( L\_{\max} \\)、\\( N\_{\mathrm{rad}} \\)、同じ波数格子のとき）。

## 〔未確定〕

- 1.3 節の \\( V^{a}\_H \\) の数値積分と、1.4 節・3.5 節の球 Bessel 変換の求積の細部（既存実装はそれぞれ専用の求積を使う）。
  結果を許容誤差の範囲で変えるだけと考えていますが、fixture との照合で確かめるまでは確定しません。
- 4 節の \\( \tilde P^{ab}\_L \\) を、選んだ基底関数の組ごとに作るか、PAO ファイルの全動径関数について作るか（結果は同じで、効率だけの問題）。
- 軌道最適化（縮約）を使う計算で、縮約後の基底関数が展開の空間に入る条件（縮約に使う原始関数の数が \\( N\_{\mathrm{rad}} \\) 以下であること）が常に満たされるか。
- 非共線スピン・スピン軌道相互作用での扱い（\\( V^{\mathrm{na}} \\) はスピンによらないので、スピンの対角成分に同じものを加えるだけと考えていますが未確認）。
- 応力への寄与（[全エネルギーと力](energy-force.md) の応力の節に、射影子展開の項は二中心の因子ごとに扱うとあるが、打ち切り因子 \\( d \\) と 4 節の直接計算の項の扱いは未確認）。
