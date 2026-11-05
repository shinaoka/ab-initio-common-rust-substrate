# 非共線・スピン軌道計算のエネルギー・力・応力

> 状態：初版（レビュー待ち）。[全エネルギーと力](energy-force.md) が共線スピン（実数の行列）に限って書いた内容を、
> 非共線スピン（`scf.SpinPolarization NC`）とスピン軌道相互作用（`scf.SpinOrbit.Coupling on`）へ広げます。
> 式は技術ノート 2（Ver. 1.0、以下 TN2）、技術ノート 1（Ver. 1.3、以下 TN1）、Kübler ら（J. Phys. F **18**, 469 (1988)）と
> マニュアル 4.0（37・38・39・41・43 節）を、このノートの記法で導出し直したものです。
> 公開資料に書かれていない規則は既存実装の挙動として〔規約〕で補いました。
> 代数的な恒等式は `tools/check_noncollinear_energy_force.py` で数値的に確かめています。

## この章の範囲〔設計〕

- 扱うもの：スピン 2×2 の密度行列とエネルギー重み付き密度行列、各エネルギー成分の評価、力、応力の可否、
  非共線の交換相関ポテンシャルの力への寄与、スピンの向きの拘束と Zeeman 項のエネルギー、
  初期スピン配置の入力、局所スピン・軌道モーメントの出力。
- 共線の場合と同じ部分（遮蔽された原子間項、原子密度補正、格子の扱い、符号と単位）は繰り返さず、
  [全エネルギーと力](energy-force.md) を参照します。
- 非共線の格子上の密度と XC ポテンシャルの式そのものは [実空間グリッド・Poisson・交換相関](grid-poisson-xc.md) の 2.4・4.6 節、
  スピン軌道の行列要素は [非局所擬ポテンシャル](nonlocal-pseudopotential.md) の 5 節、
  スピノル行列の並びと Bloch 和は [k 点と Bloch 和](kpoints-bloch.md)、密度行列の保持と初期密度は
  [密度行列と電荷ミキシング](density-mixing.md) が担当します。
- DFT+U（`Uhub`）、スピンの向きの拘束、Zeeman 項、軌道モーメント、初期磁気配置の詳細は [DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) が担当します。
  ここではエネルギー・力の縮約の観点から必要な範囲だけを要約し、式の記号をそろえます。

## 記号〔設計〕

[全エネルギーと力](energy-force.md) の記号に加えて次を使います。

| 記号 | 意味 |
|---|---|
| \\( \sigma, \sigma' \in \\{\uparrow,\downarrow\\} \\) | 全体座標の \\( z \\) 軸を量子化軸とするスピン |
| \\( c\_{\sigma I\alpha,\nu}(\mathbf k) \\) | スピノル固有ベクトル（長さ \\( 2N \\)、前半が \\( \uparrow \\)、後半が \\( \downarrow \\)） |
| \\( \rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | スピン 2×2 の実空間密度行列（複素数） |
| \\( \varepsilon^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | 同じ形のエネルギー重み付き密度行列 |
| \\( X^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | スピノル演算子 \\( \hat X \\) の行列要素 \\( \langle I\alpha\sigma,\mathbf 0\rvert\hat X\lvert J\beta\sigma',\mathbf L\rangle \\) |
| \\( \underline n(\mathbf r) \\) | 格子上の 2×2 スピン密度（成分 \\( n\_{\sigma\sigma'} \\)） |
| \\( \mathbf m(\mathbf r) \\)、\\( \hat{\mathbf m} \\) | スピン密度ベクトルとその向き \\( (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) \\) |
| \\( n'\_\uparrow, n'\_\downarrow \\) | 局所スピン座標系の密度 \\( \frac12(n \pm \lVert\mathbf m\rVert) \\) |
| \\( \underline V(\mathbf r) \\) | 格子上の 2×2 有効ポテンシャル |
| \\( \boldsymbol\sigma = (\sigma\_x,\sigma\_y,\sigma\_z) \\)、\\( \sigma\_0 \\) | Pauli 行列と 2×2 の単位行列 |
| \\( \langle \rho, X\rangle \\) | 1.4 節で定義する縮約 |

「向き」は、マニュアル 37 節のとおり、電子の磁気モーメントではなくスピン密度 \\( \mathbf m \\) の向きを指します
（物理的な磁気モーメントは逆向き）。

## 1. スピン 2×2 の密度行列

### 1.1 定義〔公開〕〔規約〕

〔公開〕TN2 Eq. 1–3：スピノル \\( \psi\_\nu = (\psi\_{\nu\uparrow}, \psi\_{\nu\downarrow}) \\) の占有和から密度演算子を作ります。
このノートの位相の規約（[疎ブロック形式](../design/sparse-block-format.md)、[k 点と Bloch 和](kpoints-bloch.md)）で、実空間の密度行列を

$$
\rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu(\mathbf k)\\, e^{+i\mathbf k\cdot\mathbf L}\\,
\overline{c\_{\sigma I\alpha,\nu}(\mathbf k)}\\, c\_{\sigma' J\beta,\nu}(\mathbf k)
$$

$$
\varepsilon^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu(\mathbf k)\\, e\_\nu(\mathbf k)\\, e^{+i\mathbf k\cdot\mathbf L}\\,
\overline{c\_{\sigma I\alpha,\nu}(\mathbf k)}\\, c\_{\sigma' J\beta,\nu}(\mathbf k)
$$

と定めます。\\( f\_\nu \le 1 \\)（スピノル 1 本に電子 1 個）で、共線の \\( g\_s \\) に当たる係数は 1 です。

〔規約〕既存実装の密度行列はこの式そのもの（\\( \bar c\\, c \\) の順、位相 \\( e^{+i\mathbf k\cdot\mathbf L} \\)、実部と虚部を別々に保持）です。
Hermite 性

$$
\rho^{\sigma'\sigma}\_{J\beta,I\alpha}(-\mathbf L) = \overline{\rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L)}
$$

から \\( \rho^{\downarrow\uparrow} \\) は \\( \rho^{\uparrow\downarrow} \\) で決まり、保持する成分は
[密度行列と電荷ミキシング](density-mixing.md) の「非共線の場合」のとおり、\\( \mathrm{Re}\\,\rho^{\uparrow\uparrow} \\)、\\( \mathrm{Re}\\,\rho^{\downarrow\downarrow} \\)、
\\( \mathrm{Re}\\,\rho^{\uparrow\downarrow} \\)、\\( \mathrm{Im}\\,\rho^{\uparrow\downarrow} \\) と、時間反転を破る項があるときの \\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} \\)、\\( \mathrm{Im}\\,\rho^{\downarrow\downarrow} \\) です。

### 1.2 時間反転と虚部〔公開〕（導出）

基底関数は実数なので、時間反転演算子は \\( \Theta = (i\sigma\_y)K \\)（\\( K \\) は複素共役）として軌道の添字を変えずに働きます。
Hamiltonian が時間反転対称で、Kramers 対がそろって占有されている（非磁性の）とき、\\( \Theta\rho\Theta^{-1} = \rho \\) と Hermite 性から

$$
\rho^{\uparrow\uparrow}\_{ab}(\mathbf L) = \overline{\rho^{\downarrow\downarrow}\_{ab}(\mathbf L)},\qquad
\rho^{\uparrow\downarrow}\_{ab}(\mathbf L) = -\rho^{\uparrow\downarrow}\_{ba}(-\mathbf L)
$$

が成り立ちます（\\( a = I\alpha \\)、\\( b = J\beta \\)）。つまり非磁性でも、スピン軌道相互作用があれば
\\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} = -\mathrm{Im}\\,\rho^{\downarrow\downarrow} \ne 0 \\) や、軌道について反対称な \\( \rho^{\uparrow\downarrow} \\) が現れます。
スピン軌道相互作用がなく、スピン依存の項がすべて軌道について実対称な行列とスピン行列の積なら、
\\( \rho^{\sigma\sigma'}\_{ab}(\mathbf L) = \rho^{\sigma\sigma'}\_{ba}(-\mathbf L) \\)（軌道について対称）です。
どちらも `tools/check_noncollinear_energy_force.py` の検査 5・6 で確かめています。

### 1.3 格子上のスピン密度との対応〔公開〕〔規約〕

〔公開〕TN2 Eq. 3、33–35：格子上の 2×2 密度は \\( n\_{\sigma\sigma'}(\mathbf r) = \sum\_\nu f\_\nu\\, \psi\_{\nu\sigma}(\mathbf r)\\,\overline{\psi\_{\nu\sigma'}(\mathbf r)} \\)、
\\( \underline n = \frac12(n\\,\sigma\_0 + \mathbf m\cdot\boldsymbol\sigma) \\) です（[グリッド](grid-poisson-xc.md) 2.4 節と同じ。\\( m\_y = -2\\,\mathrm{Im}\\,n\_{\uparrow\downarrow} \\)）。

1.1 節の \\( \rho \\) で書くと、スピンの添字が入れ替わります。

$$
n\_{\sigma\sigma'}(\mathbf r) = \sum\_{\mathbf R}\sum\_{\mathbf L}\sum\_{ab} \rho^{\sigma'\sigma}\_{ab}(\mathbf L)\\,
\phi\_a(\mathbf r - \mathbf R\_I - \mathbf R)\\,\phi\_b(\mathbf r - \mathbf R\_J - \mathbf L - \mathbf R)
$$

（\\( \mathbf R \\) は周期像、和は原子対の範囲）。対角成分は \\( \rho^{\sigma\sigma} \\) の実部だけで決まり（虚部は \\( a\leftrightarrow b \\) について反対称なので消える）、
非対角成分は \\( n\_{\uparrow\downarrow} = \overline{\sum \rho^{\uparrow\downarrow}\phi\phi} \\) です。

〔規約〕既存実装は \\( \sum\rho^{\uparrow\downarrow}\phi\phi \\) を作ってから虚部の符号を反転して \\( n\_{\uparrow\downarrow} \\) とします。上の式と同じです。

### 1.4 縮約の規則〔公開〕（導出）〔規約〕

スピノル演算子 \\( \hat X \\) の行列要素 \\( X^{\sigma\sigma'}\_{ab}(\mathbf L) \\)（Bloch 和は \\( X(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L}X(\mathbf L) \\)）について、

$$
\sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu\\, c\_\nu^\dagger X(\mathbf k)\\, c\_\nu
{}= \sum\_{\sigma\sigma'}\sum\_{\mathbf L}\sum\_{ab} \rho^{\sigma\sigma'}\_{ab}(\mathbf L)\\, X^{\sigma\sigma'}\_{ab}(\mathbf L)
{}\equiv \langle\rho, X\rangle
$$

が成り立ちます（複素共役を付けない積の和であることに注意）。\\( X \\) が Hermite なら値は実数で、Hermite 性を使って \\( \downarrow\uparrow \\) を \\( \uparrow\downarrow \\) にまとめると

$$
\langle\rho, X\rangle = \sum\_{\mathbf L, ab}\Bigl[
\sum\_{\sigma}\bigl(\mathrm{Re}\\,\rho^{\sigma\sigma}\\,\mathrm{Re}\\,X^{\sigma\sigma} - \mathrm{Im}\\,\rho^{\sigma\sigma}\\,\mathrm{Im}\\,X^{\sigma\sigma}\bigr)
{}+ 2\bigl(\mathrm{Re}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Re}\\,X^{\uparrow\downarrow} - \mathrm{Im}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Im}\\,X^{\uparrow\downarrow}\bigr)\Bigr]
$$

です（和は保持しているすべての原子対と格子ベクトル）。スピンに依存しない実対称な \\( X = O\otimes\sigma\_0 \\)（\\( T \\)、\\( S \\)、\\( V^{\mathrm{na}} \\)、スカラーの \\( V^{\mathrm{NL}} \\)）では

$$
\langle\rho, O\otimes\sigma\_0\rangle = \sum\_{\mathbf L,ab}\mathrm{Re}\bigl(\rho^{\uparrow\uparrow}\_{ab} + \rho^{\downarrow\downarrow}\_{ab}\bigr)(\mathbf L)\\, O\_{ab}(\mathbf L)
$$

に簡約されます（\\( \mathrm{Im}\\,\rho^{\sigma\sigma} \\) の寄与は Hermite 性で消える）。

〔規約〕既存実装のエネルギーと力の縮約は、すべてこの 2 つの形のどちらかです（`Re·Re − Im·Im`、非対角スピンに係数 2）。
数値検査：検査 1（バンドエネルギーとの一致、分解形の一致）。

## 2. エネルギー成分

### 2.1 一覧〔公開〕〔規約〕

出力のラベルと順序は共線の場合と同じです（[全エネルギーと力](energy-force.md) の「出力されるエネルギー成分」）。非共線での定義は次のとおりです。

| ラベル | 非共線での定義 | 備考 |
|---|---|---|
| `Uele` | \\( \sum\_{\mathbf k}w\_{\mathbf k}\sum\_\nu f\_\nu e\_\nu \\) | スピノルの固有値の和（係数 1） |
| `Ukin` | \\( \langle\rho, T\otimes\sigma\_0\rangle \\) | \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) だけ |
| `Una` | \\( \langle\rho, V^{\mathrm{na}}\otimes\sigma\_0\rangle \\)、または \\( \Delta V\sum\_p n\\, V^{\mathrm{na}} \\) | \\( n = n\_{\uparrow\uparrow}+n\_{\downarrow\downarrow} \\) |
| `Unl` | \\( \langle\rho, V^{\mathrm{NL}}\rangle \\) | スピン軌道ありでは 4 ブロックすべて（2.2 節） |
| `UH0`、`Ucore` | 共線と同じ | 原子位置だけの関数 |
| `UH1` | \\( \frac12\Delta V\sum\_p (n - n^{a})\\,\delta V\_H \\) | \\( n \\) は全密度 |
| `Uxc0`、`Uxc1` | 局所スピン座標系の多数・少数スピンの分（2.3 節） | 全体座標の \\( \uparrow\downarrow \\) ではない |
| `Uhub` | DFT+U | [DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) |
| `Ucs` | スピンの向きの拘束（2.4 節） | `scf.Constraint.NC.Spin` |
| `Uzs`、`Uzo` | スピン・軌道モーメントの Zeeman 項（2.4 節） | `scf.NC.Zeeman.Spin`、`scf.NC.Zeeman.Orbital` |
| `Uef`、`UvdW`、`Uch` | 共線と同じ（電場・分散力・内殻正孔） | `Uch` の行列は 1.4 節の一般の縮約 |

`Utot` は共線と同じ和です。

### 2.2 非局所項とスピン軌道相互作用〔公開〕〔規約〕

〔公開〕j 依存の擬ポテンシャル（TN2 Eq. 40–58）の行列要素 \\( V^{\mathrm{NL},\sigma\sigma'}\_{ab}(\mathbf L) \\) は
[非局所擬ポテンシャル](nonlocal-pseudopotential.md) の 5.3 節の式で作ります。対角スピンブロックの虚部は \\( L\_z \\) の項、
\\( \uparrow\downarrow \\) ブロックは \\( L\_- \\) の項です。

〔規約〕

- `Unl` は 1.4 節の一般形 \\( \langle\rho, V^{\mathrm{NL}}\rangle \\) で、スピン軌道の部分（対角スピンの虚部と \\( \uparrow\downarrow \\) ブロック）を含みます。
- `Unl` に入るのは擬ポテンシャルの非局所部分だけです。DFT+U、拘束、Zeeman 項のポテンシャルは、
  Hamiltonian では同じ虚部の行列に足されますが、エネルギーでは `Uhub`、`Ucs`、`Uzs`、`Uzo` として別に評価します。
- スピン軌道相互作用なし（`scf.SpinOrbit.Coupling off`、または j 依存のファイルがない）の非共線計算では、
  \\( V^{\mathrm{NL}} \\) はスカラー（j 平均したもの、非局所の章 4 節）で、\\( \langle\rho, V^{\mathrm{NL}}\otimes\sigma\_0\rangle \\) です。

### 2.3 交換相関：局所スピン座標系〔公開〕〔規約〕

〔公開〕Kübler ら（1988）と TN2 Eq. 4–5、17–25：交換相関エネルギーは、各点で 2×2 の密度を対角化した
局所座標系の密度 \\( n'\_\uparrow, n'\_\downarrow \\) に共線の汎関数を当てはめたものです。部分内殻補正を含めて

$$
E\_{\mathrm{xc}} = \int d^3r\\; \bigl(n + n^{\mathrm{pcc}}\bigr)\\,
\epsilon\_{\mathrm{xc}}\bigl(n'\_\uparrow + n^{\mathrm{pcc}}\_\uparrow,\\; n'\_\downarrow + n^{\mathrm{pcc}}\_\downarrow;\\; \nabla n'\_\uparrow, \nabla n'\_\downarrow\bigr)
$$

（GGA の勾配の扱いは [グリッド](grid-poisson-xc.md) 4.6 節：\\( n'\_\sigma \\) を格子上のスカラー場とみなす）。
\\( n'\_\sigma \\) は \\( \underline n \\) の固有値なので、全体のスピン回転で不変です。

〔規約〕

- `Uxc0`、`Uxc1` は、共線の式（[全エネルギーと力](energy-force.md) の「格子上で評価する項」）の \\( n\_\uparrow, n\_\downarrow \\) を
  局所座標系の \\( n'\_\uparrow, n'\_\downarrow \\)（\\( n'\_\uparrow \ge n'\_\downarrow \\)）に置き換えたものです：
  \\( U\_{\mathrm{xc}0} = \Delta V\sum\_p (n'\_\uparrow + n^{\mathrm{pcc}}\_\uparrow)\\,\epsilon\_{\mathrm{xc}} \\)（＋原子密度補正の半分）、`Uxc1` は \\( n'\_\downarrow \\) で同様。
- \\( n^{\mathrm{pcc}}\_\uparrow = n^{\mathrm{pcc}}\_\downarrow = \frac12 n^{\mathrm{pcc}} \\) です。`<scf.pcc.opencore` で \\( +1 \\)（\\( -1 \\)）の種は、
  全部を局所座標系の多数（少数）スピンに加えます。
- 原子密度補正 \\( \Delta E^{(0)} \\) はスピンによらないので、共線と同じく半分ずつ加えます。
- したがって `Uxc0`、`Uxc1` を個別に共線計算と比べられるのは、すべての点で \\( \mathbf m \\) が \\( +z \\) 向き（\\( n\_{\uparrow\uparrow} \ge n\_{\downarrow\downarrow} \\)）のときだけです。
  和 `Uxc0 + Uxc1` はいつでも比べられます。
- \\( \lVert\mathbf m\rVert < 10^{-14} \\) の点の扱いは [グリッド](grid-poisson-xc.md) 2.4 節のとおりです。

### 2.4 拘束と Zeeman 項〔公開〕〔規約〕

〔公開〕マニュアル 41 節：`scf.Constraint.NC.Spin`（`on`、`on2`、`off`、既定 `off`）は、スピンの向きが初期の向きから
ずれるとエネルギーが上がる罰則項を加えます。`on` は向きだけ、`on2` は向きと大きさを拘束します。
強さは `scf.Constraint.NC.Spin.v`（eV、既定 0）、対象の原子は `Atoms.SpeciesAndCoordinates` の第 12 列（1 が対象）です。
マニュアル 43 節と TN2 Eq. 139–150：`scf.NC.Zeeman.Spin`・`scf.NC.Mag.Field.Spin`（Tesla）でスピンモーメントの Zeeman 項、
`scf.NC.Zeeman.Orbital`・`scf.NC.Mag.Field.Orbital` で軌道モーメントの Zeeman 項を加えます。
磁場の向きは原子ごとの第 8・9 列（スピン）、第 10・11 列（軌道）で、第 12 列が 1 の原子にだけ働きます。

〔規約〕詳細は [DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) の 5〜7 節にあり、ここでは要約します。これらの項は、原子ごとのスピン 2×2 占有 \\( \underline N\_I \\)（格子上の密度と同じ向きの規約。同章の \\( T^{(I)} \\)）で定義されます。

$$
A^{\sigma\sigma'}\_I = \sum\_{i\in I}\sum\_{J,\mathbf L}\sum\_{j} \rho^{\sigma\sigma'}\_{Ii,Jj}(\mathbf L)\\, S\_{Ii,Jj}(\mathbf L),\qquad
\underline N\_I = \tfrac12\bigl(\overline{\underline A\_I} + \underline A\_I^{\mathsf T}\bigr)
$$

これは DFT+U の 'dual' 占有行列を軌道についてトレースしたもので、Hermite な 2×2 行列です
（'dual' 占有の定義：Han, Ozaki, Yu, PRB **74**, 045110 (2006)）。\\( \underline N\_I \\) の固有値を \\( N^\uparrow\_I \ge N^\downarrow\_I \\)、
\\( \mathbf s\_I = (\mathrm{Re}\\,N^{\uparrow\downarrow}\_I,\ -\mathrm{Im}\\,N^{\uparrow\downarrow}\_I,\ \frac12(N^{\uparrow\uparrow}\_I - N^{\downarrow\downarrow}\_I)) \\)（\\( \hbar = 1 \\) のスピン）とします。
\\( \hat{\mathbf n}^0\_I \\) は入力の角度（5 節の反転規則を適用した後）の単位ベクトルです。

- 拘束（`on`）：

$$
U\_{\mathrm{cs}} = v\sum\_{I}\bigl\lVert \underline N\_I - \underline N^0\_I\bigr\rVert\_F^2,\qquad
\underline N^0\_I = \tfrac12\Bigl[(N^\uparrow\_I + N^\downarrow\_I)\\,\sigma\_0 + (N^\uparrow\_I - N^\downarrow\_I)\\,\hat{\mathbf n}^0\_I\cdot\boldsymbol\sigma\Bigr]
$$

  （\\( \lVert\cdot\rVert\_F \\) は Frobenius ノルム、\\( v \\) は `scf.Constraint.NC.Spin.v` を Hartree にしたもの）。
  目標 \\( \underline N^0\_I \\) は現在の大きさを保って向きだけを入力の向きにしたものです。
  ポテンシャルは \\( \underline N^0\_I \\) の \\( \underline N\_I \\) への依存（固有値を通じた依存）も微分に含めます。
- 拘束（`on2`）：\\( N^\uparrow\_I - N^\downarrow\_I \\) を入力の \\( \lvert N^{\uparrow}\_{I,\mathrm{in}} - N^{\downarrow}\_{I,\mathrm{in}}\rvert \\)（第 6・7 列）に置き換えます（全数 \\( N^\uparrow\_I + N^\downarrow\_I \\) は現在の値）。
- スピンの Zeeman 項：\\( U\_{\mathrm{zs}} = -B\_s\sum\_I \hat{\mathbf n}^0\_I\cdot\mathbf s\_I \\)。\\( B\_s \\) は `scf.NC.Mag.Field.Spin` を
  \\( 1\ \mathrm{T} = 1/(2.35051742\times10^{5}) \\) a.u. で換算したものです。スピンが \\( \hat{\mathbf n}^0\_I \\) にそろうとエネルギーが下がる符号です。
- 軌道の Zeeman 項：\\( U\_{\mathrm{zo}} = -\frac12 B\_o\sum\_I \hat{\mathbf n}^{0}\_{o,I}\cdot\mathbf M\_I \\)。\\( \mathbf M\_I \\) は 6.2 節の軌道モーメントのベクトル（同章の \\( \mathbf L\_I \\)）、
  \\( B\_o \\) は `scf.NC.Mag.Field.Orbital` を同じ換算で a.u. にしたもの、\\( \hat{\mathbf n}^{0}\_{o,I} \\) は第 10・11 列の向きです。
- 入力の条件：これらはすべて非共線計算でだけ使え、`scf.Hubbard.Occupation` が `dual`（既定）であることが必要です。
  拘束とスピンの Zeeman 項は同時に使えません。軌道の Zeeman 項はスピン軌道相互作用が必要です。
  条件を満たさない入力はエラーで停止します。

〔公開〕TN2 Eq. 148–150 のスピン成分 \\( s\_y = \frac{i}{2}(N\_{\alpha\beta} - N\_{\beta\alpha}) = -\mathrm{Im}\\,N\_{\alpha\beta} \\) は上の \\( \mathbf s\_I \\) と一致します。

### 2.5 バンドエネルギーの恒等式〔公開〕（導出）

1.4 節から、収束していなくても同じ固有ベクトルについて

$$
U\_{\mathrm{ele}} = \bigl\langle \rho,\ T\otimes\sigma\_0 + V^{\mathrm{na}}\otimes\sigma\_0 + V^{\mathrm{NL}} + H^{\mathrm{grid}} + V^{\mathrm{U}} \bigr\rangle,\qquad
H^{\mathrm{grid},\sigma\sigma'}\_{ab}(\mathbf L) = \Delta V\sum\_p \phi\_a\\, V\_{\sigma\sigma'}(\mathbf r\_p)\\,\phi\_b
$$

が成り立ちます（\\( V^{\mathrm{U}} \\) は DFT+U・拘束・Zeeman の行列、\\( \underline V = (\delta V\_H\\,[+V^{\mathrm{na}}])\sigma\_0 + \underline V\_{\mathrm{xc}} \\)）。
格子の項は 1.3 節の添字の入れ替えにより \\( \langle\rho, H^{\mathrm{grid}}\rangle = \Delta V\sum\_p \mathrm{Tr}[\underline V\\,\underline n] \\) です（検査 3）。
TN2 Eq. 6 の二重計数の形 \\( E = E\_{\mathrm{band}} - \frac12\int n V\_H - \int\mathrm{Tr}(\underline V\_{\mathrm{xc}}\underline n) + E\_{\mathrm{xc}} \\) はこの恒等式の言い換えです。

## 3. 力〔公開〕〔規約〕

### 3.1 全体の形〔公開〕（導出）

\\( \mathbf F\_K = -\partial U\_{\mathrm{tot}}/\partial\mathbf R\_K \\) です（符号・単位・出力は共線と同じ）。スピノルの一般化固有値問題
\\( H^{(2)}c = e\\,S^{(2)}c \\)、\\( S^{(2)} = \sigma\_0\otimes S \\) について、占有を固定したとき

$$
\partial\sum\_\nu f\_\nu e\_\nu = \langle\rho, \partial H\rangle - \langle\varepsilon, \partial S^{(2)}\rangle,\qquad
\langle\varepsilon, \partial S^{(2)}\rangle = \sum\_{\mathbf L,ab}\mathrm{Re}\bigl(\varepsilon^{\uparrow\uparrow}\_{ab} + \varepsilon^{\downarrow\downarrow}\_{ab}\bigr)(\mathbf L)\\,\partial S\_{ab}(\mathbf L)
$$

です（TN1 Eq. 71–75 のスピノル版。検査 2）。したがって Pulay の項に必要なのは \\( \mathrm{Re}(\varepsilon^{\uparrow\uparrow}+\varepsilon^{\downarrow\downarrow}) \\) だけです。
〔規約〕既存実装も Pulay の項にはこの和だけを使います。

力の分割は共線と同じで、各項の「スピンについての和 \\( g\_s\sum\_\sigma \\)」を 1.4 節の縮約に置き換えます。

| 項 | \\( \partial U\_{\mathrm{tot}}/\partial\mathbf R\_K \\) への寄与（非共線での形） |
|---|---|
| Pulay | \\( -\langle\varepsilon, \partial S\otimes\sigma\_0\rangle \\) |
| 運動、射影子展開の \\( V^{\mathrm{na}} \\) | \\( \langle\rho, \partial O\otimes\sigma\_0\rangle \\)、\\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) だけ |
| 非局所（スピン軌道を含む） | \\( \langle\rho, \partial V^{\mathrm{NL}}\rangle \\)（3.2 節） |
| 格子：基底関数の移動 | \\( \langle\rho, \partial H^{\mathrm{grid}}\rangle \\)（3.3 節） |
| 格子：\\( n^{a} \\)、\\( n^{\mathrm{pcc}} \\)、格子モードの \\( V^{\mathrm{na}} \\)、\\( \Delta E^{(0)} \\) | 共線と同じ（\\( n^{\mathrm{pcc}} \\) は 3.4 節） |
| `Ucore + UH0` | 共線と同じ |
| `Ucs`、`Uzs`、`Uzo`、`Uhub` | 3.5 節 |

### 3.2 非局所項とスピン軌道の力〔公開〕（導出）〔規約〕

非局所の行列要素の微分 \\( \partial V^{\mathrm{NL},\sigma\sigma'}\_{ab}/\partial\mathbf R\_X \\) は [非局所擬ポテンシャル](nonlocal-pseudopotential.md) の 6 節の式
（減衰関数の微分を含む。\\( P \\) は j 射影の 2×2×軌道の行列）で、4 つのスピンブロックごとに複素数として作ります。力はこれと \\( \rho \\) の縮約

$$
\frac{\partial U\_{\mathrm{nl}}}{\partial\mathbf R\_X} = \sum\_{\mathbf L,ab}\Bigl[\sum\_\sigma\bigl(\mathrm{Re}\\,\rho^{\sigma\sigma}\\,\mathrm{Re}\\,\partial V^{\sigma\sigma} - \mathrm{Im}\\,\rho^{\sigma\sigma}\\,\mathrm{Im}\\,\partial V^{\sigma\sigma}\bigr)
{}+ 2\bigl(\mathrm{Re}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Re}\\,\partial V^{\uparrow\downarrow} - \mathrm{Im}\\,\rho^{\uparrow\downarrow}\\,\mathrm{Im}\\,\partial V^{\uparrow\downarrow}\bigr)\Bigr]
$$

です。スピン軌道の部分の力は、\\( \mathrm{Im}\\,\rho^{\sigma\sigma} \\) と \\( \rho^{\uparrow\downarrow} \\) を通じて入ります。

〔規約〕既存実装は、スピン軌道相互作用・DFT+U・拘束・Zeeman 項のどれかがあるときにこの一般形を使い、
どれもないときは \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) だけの共線と同じ形を使います（後者では一般形と数学的に同じ）。

### 3.3 格子項と非共線 XC ポテンシャル〔公開〕（導出）〔規約〕

基底関数の移動による格子の項は、\\( \rho \\) と \\( \underline V \\) を固定した \\( \Delta V\sum\_p\mathrm{Tr}[\underline V\\,\underline n] \\) の微分です。
行と列の寄与は Hermite 性で互いに複素共役になるので

$$
\frac{\partial E}{\partial\mathbf R\_K}\Big|\_{\mathrm{basis}}
{}= -2\\,\mathrm{Re}\sum\_{\alpha\in K}\sum\_{J,\mathbf L,\beta}\sum\_{\sigma\sigma'} \rho^{\sigma\sigma'}\_{K\alpha,J\beta}(\mathbf L)\\;
\Delta V\sum\_p \nabla\phi\_{K\alpha}(\mathbf r\_p - \mathbf R\_K)\\, V\_{\sigma\sigma'}(\mathbf r\_p)\\, \phi\_{J\beta}(\mathbf r\_p - \mathbf R\_J - \mathbf L)
$$

です（共線の TN1 Eq. 67、79 の係数 2 に当たる）。\\( V\_{\sigma\sigma'} = \delta\_{\sigma\sigma'}(\delta V\_H\\,[+V^{\mathrm{na}}]\\,[+V\_{\mathrm{ef}}]) + V\_{\mathrm{xc},\sigma\sigma'} \\)、
\\( \underline V\_{\mathrm{xc}} = \bar v\\,\sigma\_0 + \Delta v\\,\hat{\mathbf m}\cdot\boldsymbol\sigma \\)（[グリッド](grid-poisson-xc.md) 4.6 節）です。

**局所座標系の角度は微分しなくてよい理由。** \\( E\_{\mathrm{xc}} \\) は \\( \underline n \\) に固有値 \\( n'\_\pm \\) を通じてだけ依存します。
2×2 Hermite 行列の固有値の微分は \\( \partial n'\_\pm/\partial n\_{\sigma'\sigma} = (P\_\pm)\_{\sigma\sigma'} \\)、\\( P\_\pm = \frac12(\sigma\_0 \pm \hat{\mathbf m}\cdot\boldsymbol\sigma) \\)（固有射影）なので

$$
\frac{\partial E\_{\mathrm{xc}}}{\partial n\_{\sigma'\sigma}(\mathbf r)} = v'\_\uparrow (P\_+)\_{\sigma\sigma'} + v'\_\downarrow (P\_-)\_{\sigma\sigma'} = V\_{\mathrm{xc},\sigma\sigma'}(\mathbf r)
$$

で、TN2 Eq. 18–19 の「局所座標系で求めて全体座標系に戻す」ポテンシャルはそのまま \\( E\_{\mathrm{xc}} \\) の厳密な微分です（検査 4）。
角度 \\( (\theta,\varphi) \\) の変化は \\( P\_\pm \\) に含まれており、力に別の項は出ません。GGA でも、\\( v'\_\sigma \\) が \\( n'\_\sigma \\) についての（離散化した）汎関数微分であれば同じです。
ただし \\( \lVert\mathbf m\rVert \to 0 \\) の点では \\( n'\_\pm \\) が \\( \underline n \\) の滑らかな関数でなく（\\( n'\_+ - n'\_- = \lVert\mathbf m\rVert \\)）、GGA の \\( \nabla n'\_\sigma \\) がそこで折れ目を持つため、
有限差分との一致がその近傍で悪くなることがあります（〔未確定〕の項目）。

〔規約〕既存実装は、格子の各点で \\( (V\_{\uparrow\uparrow}, V\_{\downarrow\downarrow}, \mathrm{Re}\\,V\_{\downarrow\uparrow}, \mathrm{Im}\\,V\_{\downarrow\uparrow}) \\) を用意し、
基底関数の移動による密度成分の変化と \\( \sum\partial n\_{\uparrow\uparrow}V\_{\uparrow\uparrow} + \partial n\_{\downarrow\downarrow}V\_{\downarrow\downarrow} + 2\\,\mathrm{Re}(\partial n\_{\uparrow\downarrow}V\_{\downarrow\uparrow}) \\) の形で縮約します。上の式と同じです。

### 3.4 部分内殻補正の力〔公開〕（導出）〔規約〕

\\( n^{\mathrm{pcc}} \\) は局所座標系の両スピンに半分ずつ入るので、厳密な微分は

$$
\frac{\partial E\_{\mathrm{xc}}}{\partial\mathbf R\_K}\Big|\_{\mathrm{pcc}} = \Delta V\sum\_p \tfrac12\bigl(v'\_\uparrow + v'\_\downarrow\bigr)\\,\frac{\partial n^{\mathrm{pcc}}\_K}{\partial\mathbf R\_K}
{}= \Delta V\sum\_p \tfrac12\mathrm{Tr}\\,\underline V\_{\mathrm{xc}}\\,\frac{\partial n^{\mathrm{pcc}}\_K}{\partial\mathbf R\_K}
$$

です。開殻内殻の種（`<scf.pcc.opencore` が \\( \pm1 \\)）では、局所座標系の多数（少数）スピンのポテンシャル \\( v'\_\uparrow \\)（\\( v'\_\downarrow \\)）が掛かるのが厳密です。

〔規約〕既存実装は通常の種で \\( \mathrm{Tr}\\,\underline V\_{\mathrm{xc}} = V\_{\mathrm{xc},\uparrow\uparrow} + V\_{\mathrm{xc},\downarrow\downarrow} \\) を使います（厳密な式と同じ）。
開殻内殻の種では、非共線でも全体座標の \\( V\_{\mathrm{xc},\uparrow\uparrow} \\)（\\( +1 \\)）または \\( V\_{\mathrm{xc},\downarrow\downarrow} \\)（\\( -1 \\)）の 2 倍に \\( \frac12\partial n^{\mathrm{pcc}}\_K \\) を掛けます。
これは \\( \hat{\mathbf m} \\) が \\( +z \\) 向きの点でだけ厳密な微分と一致します。

### 3.5 DFT+U・拘束・Zeeman 項の力〔規約〕

これらのエネルギーは \\( \rho \\) と \\( S \\) から作る占有 \\( \underline N\_I \\)（とその軌道成分）の関数です。
\\( \rho \\) の変化による分は Hamiltonian に入れたポテンシャルを通じて Pulay の項にまとまり、残りは \\( \rho \\) を固定した \\( S \\) の微分による
\\( \sum\_I\mathrm{Tr}[(\partial U/\partial\underline N\_I)\\,\partial\underline N\_I/\partial\mathbf R] \\) です。
既存実装は、DFT+U（'dual'）の力と同じ式に、拘束・Zeeman 項のポテンシャルを足した有効ポテンシャル行列を使ってこれを評価します。
式は [DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) の 11 節です（\\( \tfrac12(V\partial S + \partial S V) \\) と \\( \rho \\) の 1.4 節の縮約）。

### 3.6 既存実装の力の検査結果〔規約〕

配布物の力の検査（`openmx -forcetest`、マニュアル 72 節）の結果ファイルには、非共線の入力について解析的な力と
\\( \pm 3\times10^{-4} \\) bohr の中心差分の差が記録されています。目安として次の程度です（全項、原子 1 の各成分の最大）。

| 入力 | 内容 | 解析と数値の差（Hartree/bohr） |
|---|---|---|
| `CH_NC.dat` | 非共線、スピン軌道なし（分子） | \\( 4\times10^{-4} \\) |
| `Cr2.dat` | 非共線、スピン軌道、DFT+U、拘束（Cluster） | \\( 1\times10^{-4} \\) |
| `GaAs_NC_SO.dat` | 非共線、スピン軌道（Band） | \\( 1\times10^{-3} \\) |
| `GaAs_NC_SO+U.dat` | 非共線、スピン軌道、DFT+U | \\( 3\times10^{-4} \\) |

格子の離散化（egg-box 効果）と SCF の収束の誤差を含む値です。再実装の許容誤差の上限の参考にします。

## 4. 応力〔規約〕〔設計〕

〔規約〕既存実装は、非共線計算（`scf.SpinPolarization NC`）と応力の計算を併用できません。応力が要求されると
（`scf.stress.tensor on`、格子を動かす最適化、NPT の分子動力学）、入力の検査でエラーとして停止します。
スピン軌道相互作用の有無にはよりません。応力と併用できない組み合わせは、ほかに軌道最適化と DFT+U の 'full' 射影子です
（[全エネルギーと力](energy-force.md) の「出力と格子ベクトルの勾配」）。

〔設計〕再実装の初版も同じく非共線と応力を併用不可とし、同じ入力をエラーにします。将来対応する場合の見通しを記録しておきます。

- 行列要素の項（運動、非局所、射影子展開の \\( V^{\mathrm{na}} \\)）と Pulay の項は、共線の表（TN23 Eq. 16–31、62–67）の
  \\( g\_s\sum\_\sigma\rho^\sigma \\) を 1.4 節の縮約に置き換えるだけです。スピン軌道の項も、二中心の因子ごとの
  \\( \partial/\partial\tau\_\gamma\\,\tau\_\eta \\) で同じように扱えます。
- 交換相関の体積項と LDA の項は \\( n'\_\sigma \\) で書けば共線と同じ形です。GGA の追加分も \\( n'\_\sigma \\) の差分勾配で書けますが、
  \\( \lVert\mathbf m\rVert\to 0 \\) の点の扱いを決める必要があります。
- DFT+U・拘束・Zeeman の項の応力は未検討です。

## 5. 初期スピン配置〔公開〕〔規約〕

〔公開〕マニュアル 37 節：非共線計算では `Atoms.SpeciesAndCoordinates` の各行が 13 列です。

| 列 | 内容 |
|---|---|
| 1–2 | 通し番号、原子種 |
| 3–5 | 座標 |
| 6–7 | 初期の上向き・下向きの電子数 \\( N^{\uparrow}\_{I,\mathrm{in}}, N^{\downarrow}\_{I,\mathrm{in}} \\) |
| 8–9 | スピンの向き \\( (\theta, \varphi) \\)（度）。初期スピン密度の向き、拘束・スピン Zeeman 項の目標の向き |
| 10–11 | 軌道モーメントの向き（度）。軌道 Zeeman 項の目標の向き |
| 12 | 拘束・Zeeman 項を適用するか（1／0） |
| 13 | DFT+U の軌道偏極を促すか（`on`／`off`） |

〔規約〕

- 角度は度で読み、ラジアンに直します。\\( \theta \\)（スピン・軌道とも）の絶対値が \\( 10^{-10} \\) 度より小さいときは、
  変換前に幅 \\( 10^{-5} \\) 度の一様乱数（中心 0）を加えます。スピンの \\( \theta, \varphi \\) には、変換後さらに幅 \\( 10^{-10} \\) ラジアンの一様乱数を加えます。
  正確な対称配置からわずかにずらすためで、結果への影響は無視できます。
- \\( N^{\uparrow}\_{I,\mathrm{in}} < N^{\downarrow}\_{I,\mathrm{in}} \\) のときは、スピンの向きを反対向き \\( -\hat{\mathbf n} \\) に置き換え、
  \\( N^{\uparrow}\_{I,\mathrm{in}} \\) と \\( N^{\downarrow}\_{I,\mathrm{in}} \\) を入れ替えます。軌道モーメントの向きも同時に反対向きにします。
  新しい角度は \\( -\hat{\mathbf n} \\) の極座標（\\( \theta\in[0,\pi] \\)、\\( \varphi\in[-\pi/2, 3\pi/2] \\)、6.3 節）です。
  例：第 6・7 列が 6.0、8.0、角度が \\( (42^\circ, 79^\circ) \\) の原子は、\\( N^\uparrow = 8 \\)、\\( N^\downarrow = 6 \\)、向き \\( (138^\circ, 259^\circ) \\) として扱われます。
  これにより、以後は常に「上向き」が入力の向きに沿った多数スピンになります。拘束・Zeeman 項の目標の向きも反転後の向きです。
- 初期密度は、原子 \\( I \\) の原子密度を \\( (N^\uparrow\_{I,\mathrm{in}} + N^\downarrow\_{I,\mathrm{in}})/Z\_I \\) 倍した全密度と、
  \\( (N^\uparrow\_{I,\mathrm{in}} - N^\downarrow\_{I,\mathrm{in}})/Z\_I \\) 倍した大きさを入力の向きに持つ \\( \mathbf m \\) から作ります
  （[密度行列と電荷ミキシング](density-mixing.md) の 6 節）。
- 第 8 列以降が欠けた行の扱いと軌道分極フラグは [DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) の 8 節を参照してください。

## 6. 局所スピン・軌道モーメントの出力〔公開〕〔規約〕

### 6.1 原子ごとのスピンモーメント〔規約〕

原子ごとの 2×2 の Mulliken 型の占有 \\( \underline M\_I \\) を、2.4 節の \\( \underline A\_I \\) から次のように作ります。

$$
M^{\sigma\sigma}\_I = \mathrm{Re}\\,A^{\sigma\sigma}\_I,\qquad M^{\uparrow\downarrow}\_I = \overline{A^{\uparrow\downarrow}\_I},\qquad M^{\downarrow\uparrow}\_I = A^{\uparrow\downarrow}\_I
$$

（\\( \underline A\_I \\) の行は原子 \\( I \\) の軌道、列は相手の全軌道。軌道ごとの分解では \\( i\in I \\) の和をとらない）。\\( \underline M\_I \\) を
1.3 節と同じ規則 \\( \mathbf m = (2\\,\mathrm{Re}\\,M^{\uparrow\downarrow}, -2\\,\mathrm{Im}\\,M^{\uparrow\downarrow}, M^{\uparrow\uparrow}-M^{\downarrow\downarrow}) \\) でベクトルにし、
\\( N^\uparrow\_I, N^\downarrow\_I = \frac12(\mathrm{Tr}\\,\underline M\_I \pm \lVert\mathbf m\_I\rVert) \\) と \\( \mathbf m\_I \\) の向き \\( (\theta\_I, \varphi\_I) \\) を求めます。

- この \\( \underline M\_I \\) は 2.4 節の \\( \underline N\_I \\)（Hermite 部分）とは異なり、非対角成分に原子 \\( I \\) の行だけを使います。
  両者の差は \\( \rho^{\uparrow\downarrow} \\) の軌道について反対称な部分から来ます。1.2 節により、スピン軌道相互作用のない系では差は 0 です。
  スピン軌道相互作用のある系では、対称性の低い原子で差が生じえて、時間反転対称な系でも原子ごとに見かけのモーメントが出る可能性があります
  （単純なモデルで確認：検査 6。GaAs のように原子位置の対称性が高い場合は対称性で 0 になる）。全原子の和は両者で一致します。
- 出力のモーメントは \\( N^\uparrow\_I - N^\downarrow\_I \\)（\\( \mu\_B \\)、\\( g = 2 \\) とみなす）で、常に 0 以上です。向きは角度で表します。
- 系全体のスピンモーメントは \\( \sum\_I (N^\uparrow\_I - N^\downarrow\_I)\\,\hat{\mathbf m}\_I \\) の大きさと向きです。
- 求めた \\( N^\uparrow\_I, N^\downarrow\_I, \theta\_I, \varphi\_I \\) は、MD・構造最適化の次のステップの初期値などにも使われます
  （拘束・Zeeman 項の目標の向きは入力の値のまま変わりません）。

### 6.2 軌道モーメント〔公開〕〔規約〕

〔公開〕マニュアル 39 節、TN2 Eq. 151、160–162：非共線計算では原子ごとの軌道モーメントを、局在基底での軌道角運動量の期待値として求めます。
実基底では角運動量の行列要素は純虚数なので、スピン軌道相互作用がなければ 0 です。

〔規約〕既存実装は、原子 \\( I \\) 上の \\( \langle\phi\_{Im}\rvert\hat{\mathbf L}\lvert\phi\_{In}\rangle \\)（一中心の行列）と、
\\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} \\)、\\( \mathrm{Im}\\,\rho^{\downarrow\downarrow} \\) の 'dual' 型の積（\\( \rho S \\) と \\( S\rho \\) の平均）から \\( \mathbf M\_I \\) を作り、大きさと角度を出力します。
式と符号は [DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) の 7 節です。角度の規約は下の 6.3 節と同じです。

### 6.3 出力の書式〔公開〕〔規約〕

〔公開〕マニュアル 37・39 節：`System.Name.out` の Mulliken の節に、全スピンモーメントとその角度（度）、原子ごとの
`Up`、`Down`、`Sum`、`Diff`、`theta`、`phi` の表、和の行が出ます。続く軌道ごとの分解にも角度が付きます。
軌道モーメントの節には、全軌道モーメントと角度、原子ごとの値と角度、軌道ごとの分解が出ます。

〔規約〕

- `Up`、`Down` は \\( N^\uparrow\_I, N^\downarrow\_I \\)、`Sum` は和、`Diff` は差、`theta`・`phi` は \\( \hat{\mathbf m}\_I \\) の角度（度）です。
  先頭行の「Total spin moment (muB)」は 6.1 節の全スピンモーメントです。
- 角度は \\( \theta\in[0,\pi] \\)、\\( \varphi \\) は \\( x \ge 0 \\) で \\( \arcsin(y/\sqrt{x^2+y^2}) \\)、\\( x < 0 \\) で \\( \pi - \arcsin(\cdot) \\)、
  すなわち \\( [-\pi/2, 3\pi/2] \\) の範囲で求めます。大きさが \\( 10^{-14} \\) 未満なら \\( \theta = \pi/2 \\)、\\( \varphi = 0 \\)、\\( xy \\) 成分が \\( 10^{-14} \\) 未満なら \\( \varphi = 0 \\) です。
  度に直した値を 360 で割った余り（符号付き）として書くので、\\( \varphi \\) は \\( -90 \\) から \\( 270 \\) の範囲の値になります。
- SCF の各反復の標準出力には、原子ごとに（既定の出力量では最初の 20 原子）`MulP` の上・下・和・差と \\( (\theta,\varphi) \\)、
  軌道モーメントの大きさ `Ml` と角度、スピンと軌道を合わせたベクトル \\( (N^\uparrow\_I - N^\downarrow\_I)\hat{\mathbf m}\_I + \mathbf M\_I \\) の大きさ `Ml+s` と角度が出ます。

## 7. 設計〔設計〕

- **密度行列の型。** 再実装では \\( \rho^{\sigma\sigma'}(\mathbf L) \\)、\\( \varepsilon^{\sigma\sigma'}(\mathbf L) \\) を原子対ごとの複素数の 2×2 スピンブロック
  （[疎ブロック形式](../design/sparse-block-format.md) の非共線の形）で持ち、\\( \downarrow\uparrow \\) も保持します。
  既存実装の 4〜6 個の実数配列は、fixture の入出力のときだけ変換します。
- **縮約は 1 つの関数。** エネルギー・力・Pulay のすべてを 1.4 節の \\( \langle\rho, X\rangle \\)（\\( X \\) は 2×2 ブロックまたは \\( O\otimes\sigma\_0 \\)）で書き、
  共線の場合は \\( \rho = \rho^\sigma\oplus\rho^{\sigma'} \\) の対角スピンだけの特別な場合として同じ関数を通します。
  スピンに依存しない \\( X \\) では \\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) の近道を使ってよいが、結果は一般形と一致させます。
- **格子上の密度の向き。** 格子では \\( \underline n \\)（1.3 節の規約、\\( n\_{\uparrow\downarrow} = (m\_x - i m\_y)/2 \\)）を保持し、
  \\( \rho \\) からの組み立てで添字を入れ替えます。局所座標系への変換はエネルギー・XC の評価の中だけで行い、格子上の保持形式は変えません。
- **`Uxc0`・`Uxc1` の分け方。** fixture と比べるため、既存実装と同じく局所座標系の多数・少数スピンで分けます。
- **開殻内殻の力。** 既定は厳密な微分（3.4 節の \\( v'\_\uparrow \\) または \\( v'\_\downarrow \\)）とし、既存実装と比べるための互換モード
  （全体座標の対角成分）を用意します。
- **Mulliken のスピンモーメント。** 既定は既存実装と同じ定義（6.1 節）とし、Hermite 部分による定義を選択肢として持ちます。
- **応力。** 非共線と応力の併用は初版ではエラーとします（4 節）。

## 8. テスト〔設計〕

### 8.1 代数の検査（`tools/check_noncollinear_energy_force.py`）

| 番号 | 内容 | 結果の目安 |
|---|---|---|
| 1 | スピノルのバンドエネルギー \\( = \langle\rho, H\rangle \\)、実部・虚部の分解形、\\( \rho \\) の Hermite 性 | \\( 10^{-10} \\) |
| 2 | Pulay の恒等式（スピノル、\\( S^{(2)} = \sigma\_0\otimes S \\)）と、\\( \mathrm{Re}(\varepsilon^{\uparrow\uparrow}+\varepsilon^{\downarrow\downarrow}) \\) だけで足りること | \\( 10^{-7} \\)（差分）、\\( 10^{-12} \\)（2 つの形の差） |
| 3 | \\( \Delta V\sum\mathrm{Tr}[\underline V\\,\underline n] = \langle\rho, H^{\mathrm{grid}}\rangle \\)、\\( \underline n \\) の Hermite 性（1.3 節の添字の入れ替え） | \\( 10^{-10} \\) |
| 4 | 局所座標系で求めた \\( \underline V\_{\mathrm{xc}} \\) が \\( E\_{\mathrm{xc}} \\) の微分（モデル汎関数、Hermite な摂動の差分） | \\( 10^{-7} \\) |
| 5 | スピン軌道なしの全体スピン回転：固有値、\\( \mathrm{Re}(\rho^{\uparrow\uparrow}+\rho^{\downarrow\downarrow}) \\) の縮約が不変、原子モーメントがベクトルとして回る | \\( 10^{-10} \\) |
| 6 | 時間反転対称（スピン軌道あり）：\\( \rho^{\uparrow\downarrow} \\) が軌道について反対称、Hermite 部分の原子モーメントが 0 | \\( 10^{-8} \\) |
| 7 | 拘束の目標 \\( \underline N^0 \\) の 2 通りの書き方の一致と向きの復元 | 丸め誤差 |
| 8 | \\( H^{\uparrow\downarrow} = 0 \\) で実数の \\( H^{\sigma\sigma} \\) なら、スピノルの結果が共線 2 つと一致 | \\( 10^{-10} \\) |

### 8.2 計算全体の不変量

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 共線への帰着 | スピン軌道なし、全原子 \\( \theta = 0 \\) の `NC` と、同じ入力の `On`：`Utot` と全成分（`Uxc0`、`Uxc1` を含む。全点で \\( n\_{\uparrow\uparrow}\ge n\_{\downarrow\downarrow} \\) となる系で）、力、モーメント | \\( 10^{-8} \\) Hartree、\\( 10^{-6} \\) Hartree/bohr（SCF の収束の範囲） |
| 非磁性への帰着 | スピン軌道なし、\\( N^\uparrow = N^\downarrow \\) の `NC` と `Off` | 同上 |
| 全体スピン回転 | スピン軌道なしで、全原子の向きを同じ回転で回す：`Utot`、各成分、力が不変、モーメントの向きが回る | 同上、角度 \\( 10^{-4} \\) 度 |
| 時間反転 | スピン軌道ありで、全原子の向きを反転（\\( \theta\to\pi-\theta \\)、\\( \varphi\to\varphi+\pi \\)）：`Utot`、力が不変、スピン・軌道モーメントが反転 | 同上 |
| 回転の非不変性 | スピン軌道ありで向きを回すとエネルギーが変わること（磁気異方性。符号と大きさは fixture で確認） | — |
| 恒等式 | \\( U\_{\mathrm{ele}} = \langle\rho, H\rangle \\)（2.5 節） | \\( 10^{-10} \\) |
| 力の総和 | 行列要素の項（スピン軌道を含む）・Pulay の項の総和が 0 | \\( 10^{-10} \\) |
| 拘束 | `scf.Constraint.NC.Spin.v` を 0 にすると拘束なしと一致、目標の向きにそろった解で `Ucs` が 0 | \\( 10^{-10} \\) |

### 8.3 有限差分

- 原子 1 個を \\( \pm 3\times10^{-4} \\) bohr 動かした `Utot` の中心差分と解析的な力（全項）。
  対象：スピン軌道ありの小さな分子（例：Pt₂、`Cr2.dat` に相当する Cr₂）と、非共線の結晶（`GaAs_NC_SO.dat` に相当）。
  許容誤差は既存実装の水準（3.6 節）を上限とし、\\( 10^{-3} \\) Hartree/bohr から始めます。
- 項ごとの検査：\\( \rho \\)、\\( \varepsilon \\) を固定した行列要素の項（`Ukin`、`Unl` のスピン軌道部分を含む）と Pulay の項は \\( 10^{-7} \\)。
- 格子の項は、\\( \lVert\mathbf m\rVert \\) が 0 になる点を含まない配置（例：強磁性的にそろった配置）と含む配置の両方で調べ、
  GGA の折れ目の影響を記録します。

### 8.4 fixture との照合

既存実装から次を書き出して比べます。

- エネルギー：`Uele`、`Ukin`、`UH0`、`UH1`、`Una`、`Unl`、`Uxc0`、`Uxc1`、`Ucore`、`Uhub`、`Ucs`、`Uzs`、`Uzo`、`Utot`。
- 中間量：収束後の \\( \rho^{\sigma\sigma'}(\mathbf L) \\)（6 成分）、\\( \mathrm{Re}\\,\varepsilon^{\sigma\sigma} \\)、スピン軌道を含む \\( V^{\mathrm{NL},\sigma\sigma'}(\mathbf L) \\)、
  格子上の \\( \underline n \\)（全体座標）と \\( \underline V\_{\mathrm{xc}} \\)。これらから固有値問題を解かずにエネルギー成分と力を再計算できます。
- 出力：Mulliken の節の \\( N^\uparrow, N^\downarrow, \theta, \varphi \\)、全スピンモーメント、軌道モーメント、`<coordinates.forces`。
- 対象：配布物の入力例 `Cr2.dat`（Cluster、スピン軌道、DFT+U、拘束）と `GaAs.dat`（Band、スピン軌道）、
  力の検査用の `CH_NC.dat`、`GaAs_NC_SO.dat`、`GaAs_NC_SO+U.dat`。

| 量 | 許容誤差の目安（同じ格子・同じ基底） |
|---|---|
| 各エネルギー成分、`Utot` | \\( 10^{-6} \\) Hartree |
| \\( N^\uparrow, N^\downarrow \\)、モーメント | \\( 10^{-5} \\) |
| 角度 | \\( 10^{-3} \\) 度（モーメントが \\( 0.1\\,\mu\_B \\) 以上の原子だけ） |
| 力 | \\( 10^{-5} \\) Hartree/bohr |

## 〔未確定〕

- `Uzo` の力が、軌道モーメントの \\( S \\) を通じた依存（dual 型の積）を含むかどうか（[DFT+U・制約 DFT・磁気的な機能](dftu-magnetism.md) の 11 節の扱いとの照合）。
- 開殻内殻の種の部分内殻補正の力（3.4 節）で、既存実装の全体座標の成分を使う挙動を再実装の既定にするか。
- Mulliken のスピンモーメント（6.1 節）の行だけを使う定義が、スピン軌道ありの低対称の系で実際にどの程度の見かけのモーメントを出すか（実系での確認）。
- GGA の \\( \underline V\_{\mathrm{xc}} \\) が離散的に厳密な微分になっているか（[全エネルギーと力](energy-force.md) と共通）と、
  \\( \lVert\mathbf m\rVert\to 0 \\) の点での \\( n'\_\sigma \\) の折れ目の扱い。
- `Uxc0`・`Uxc1` の分け方で、GGA のスピンごとのエネルギー密度をどう割り振るか（共線の章と共通）。
- 角度に加える乱数の再現性（fixture を作るときに乱数の影響が許容誤差より十分小さいことの確認）。
- `DC`、`Krylov`、`NEGF`、`Cluster-LNO` などの解法での非共線のエネルギー（Harris 型を含む）と力。
- 非共線の応力（既存実装にない。4 節の見通しを実装するか）。
- （解消）[非局所擬ポテンシャル](nonlocal-pseudopotential.md) の 6 節の縮約は、この章の 1.4 節の \\( \sum\rho^{\sigma\sigma'}\_{ab}(\mathbf L)V^{\sigma\sigma'}\_{ab}(\mathbf L) \\)（複素共役なし）に統一しました。
