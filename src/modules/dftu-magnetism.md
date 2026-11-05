# DFT+U・制約 DFT・磁気的な機能

> 状態：初版（レビュー待ち）。式は技術ノート 6（LDA+U、Ver. 1.0、以下 TN6）、技術ノート 2（非共線、以下 TN2）、
> Han–Ozaki–Yu（PRB **74**, 045110 (2006)）、Ryee–Han（Sci. Rep. **8**, 9559 (2018)；JPCM **30**, 275802 (2018)）、
> マニュアル 37–43 節をこのノートの記法で書き直し、公開資料にない規約を既存実装で確認したものです。
> 式の整合性（エネルギーとポテンシャル、力、Coulomb テンソルの不変量、拘束項の閉じた形）は
> `tools/check_dftu_magnetism.py` で数値的に確かめています。

## この章の範囲

原子上の局在軌道の占有行列を通して Hamiltonian に項を加える機能をまとめます。

- DFT+U：簡略形（Dudarev 型、`scf.DFTU.Type 1`）と一般形（Liechtenstein 型、`scf.DFTU.Type 2`）、
  占有行列の 3 種類の定義（onsite・full・dual）、二重計数補正 4 種、エネルギー `Uhub`、ポテンシャル、力。
- 軌道分極の促進（`Atoms.SpeciesAndCoordinates` の軌道分極フラグ）。
- 非共線スピンの向き（と大きさ）の拘束（`scf.Constraint.NC.Spin`、エネルギー `Ucs`）。
- スピン・軌道磁気モーメントの Zeeman 項（`scf.NC.Zeeman.Spin`・`scf.NC.Zeeman.Orbital`、エネルギー `Uzs`・`Uzo`）。
- 軌道磁気モーメントの評価、第二変分法による磁気異方性エネルギー（MAE）の手順。
- 初期磁気配置の入力（座標の欄のスピン・角度・フラグ）。

格子上の非共線密度と XC は [グリッド・Poisson・交換相関](grid-poisson-xc.md)、密度行列と再開は
[密度行列と電荷ミキシング](density-mixing.md)、全エネルギーの和と力の枠組みは [全エネルギーと力](energy-force.md) にあります。
スピン軌道相互作用の行列要素そのものは扱いません。

## 記号〔設計〕

[記法](../design/notation.md)、[疎ブロック形式](../design/sparse-block-format.md)、[全エネルギーと力](energy-force.md) の記号に加えて、次を使います。

| 記号 | 意味 |
|---|---|
| \\( \rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) \\) | 実空間の密度行列。\\( \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\nu f\_\nu e^{i\mathbf k\cdot\mathbf L}\\, \overline{c\_{\sigma I\alpha,\nu}}\\, c\_{\sigma' J\beta,\nu} \\)（[密度行列](density-mixing.md)） |
| \\( s = (I, l, \zeta) \\) | シェル：原子 \\( I \\) の、角運動量 \\( l \\)・多重度 \\( \zeta \\) の \\( 2l+1 \\) 個の軌道の組 |
| \\( a, b, c, d \\) | シェル内の実球面調和関数の番号（[PAO 基底](pao-basis.md) の並び、0 から \\( 2l \\)） |
| \\( n^{\sigma\sigma'}\_{s,ab} \\) | シェル \\( s \\) の占有行列（下で定義） |
| \\( N^{\sigma\sigma'}\_s = \sum\_a n^{\sigma\sigma'}\_{s,aa} \\)、\\( N\_s = N^{\uparrow\uparrow}\_s + N^{\downarrow\downarrow}\_s \\) | シェルの占有数（スピン成分ごと、合計） |
| \\( v^{\sigma\sigma'}\_{s,ab} \\) | シェルのポテンシャル行列（\\( = \partial E / \partial n^{\sigma\sigma'}\_{s,ab} \\)） |
| \\( U\_s, J\_s \\) | シェルの U・J（入力は eV、内部は Hartree） |
| \\( U\_{abcd} = \langle ab \vert V \vert cd\rangle \\) | 実球面調和関数での Coulomb テンソル（\\( a, c \\) が電子 1、\\( b, d \\) が電子 2） |
| \\( \boldsymbol\sigma = (\sigma\_x, \sigma\_y, \sigma\_z) \\) | Pauli 行列 |
| \\( T^{(I)} \\) | 原子 \\( I \\) の 2×2 スピン密度行列（下で定義） |
| \\( \hat{\mathbf e}(\theta,\varphi) = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) \\) | 方向ベクトル |

**行列要素と微分の向き**〔設計〕。密度行列を \\( \rho\_{ab} = \sum f\\, \overline{c\_a}\\, c\_b \\) と定義しているので、バンドエネルギーは
\\( \sum\_{ab} H\_{ab}\\, \rho\_{ab} \\)（転置なしの成分ごとの積の和）です。この章では、エネルギー \\( E(n) \\) に対してポテンシャルを
\\( v\_{ab} = \partial E/\partial n\_{ab} \\)（\\( n \\) の各成分を独立とみなした偏微分）と定め、Hamiltonian に \\( v \\) をそのまま（転置せずに）加えます。
共線計算では \\( n \\) が実対称なので転置は区別できませんが、非共線計算ではこの向きを守らないと符号が反転します
（例：\\( E = \tfrac12 U\\,[\mathrm{Tr}\\,n - \mathrm{Tr}\\,n^2] \\) なら \\( v = U(\tfrac12 I - n^{\mathsf T}) \\)。TN2 の Eq. 115 の後の式と同じ）。

**実空間の行列積**〔設計〕。原子対ブロックの行列 \\( A, B \\) の積を
\\( (AB)\_{I\alpha,J\beta}(\mathbf L) = \sum\_{K,\mathbf L'}\sum\_\gamma A\_{I\alpha,K\gamma}(\mathbf L')\\, B\_{K\gamma,J\beta}(\mathbf L-\mathbf L') \\) と書きます。
これは k 空間の積 \\( A(\mathbf k)B(\mathbf k) \\) に対応します（[密度行列](density-mixing.md) の和則の節）。非共線では各スピンブロックに
スピン非依存の \\( S \\) を掛けます。\\( (\cdot)\_{II} \\) は原子 \\( I \\) の \\( \mathbf L = 0 \\) の対角ブロックを、\\( (\cdot)\_s \\) はさらにシェル \\( s \\) の部分を表します。

## 1. 占有行列〔公開〕〔規約〕

〔公開〕TN6 Eq. 5–16、TN2 Eq. 103–114：占有行列は、射影演算子の期待値

$$
n^{\sigma\sigma'}\_{s,ab} = \sum\_\nu f\_\nu \langle \psi\_\nu \vert \hat n^{\sigma\sigma'}\_{s,ab} \vert \psi\_\nu\rangle
$$

で、演算子の取り方により 3 種類あります。双対軌道 \\( \vert\tilde\phi\rangle = \sum S^{-1}\vert\phi\rangle \\) を使い、
onsite は \\( \vert\tilde\phi\_a\rangle\langle\tilde\phi\_b\vert \\)、full は \\( \vert\phi\_a\rangle\langle\phi\_b\vert \\)、
dual は \\( \tfrac12(\vert\tilde\phi\_a\rangle\langle\phi\_b\vert + \vert\phi\_a\rangle\langle\tilde\phi\_b\vert) \\) に対応します。
密度行列と重なり行列で書くと（このノートの添字の向きで）

| `scf.Hubbard.Occupation` | 占有行列 | 和則 \\( \sum\_I \mathrm{Tr}\\, n\_{II} = N\_{\mathrm e} \\) |
|---|---|---|
| `onsite` | \\( n\_{II} = \rho\_{II} \\)（\\( \mathbf L = 0 \\) の対角ブロックそのもの） | 成り立たない |
| `full` | \\( n\_{II} = (S\rho S)\_{II} \\) | 成り立たない |
| `dual`（既定） | \\( n\_{II} = \tfrac12(\rho S + S\rho)\_{II} \\) | 成り立つ |

で、シェル \\( s \\) の占有行列 \\( n\_s \\) はその \\( (2l+1)\times(2l+1) \\) の部分（非共線ではスピンの 2×2 を含めて \\( 2(2l+1) \\) 次元）です。
共線ではスピンごと（\\( \sigma = \sigma' \\)）に定義します。

〔規約〕添字の向き：既存実装の dual は \\( n^{\sigma\sigma'}\_{ab} = \tfrac12[(\rho^{\sigma\sigma'} S)\_{ab} + (S\rho^{\sigma\sigma'})\_{ab}] \\) で、TN2 Eq. 114 と同じです。
TN6 Eq. 16 はこの転置の形で書かれていますが、TN6 は共線（実対称）だけを扱うので差はありません。
非共線の dual・full の \\( n \\) は、\\( 2(2l+1) \\) 次元の行列として Hermite です（\\( n^{\downarrow\uparrow} = (n^{\uparrow\downarrow})^\dagger \\)）。
onsite の非共線では \\( n^{\uparrow\downarrow}\_{ab} = \rho^{\uparrow\downarrow}\_{Ia,Ib}(0) \\)、\\( n^{\downarrow\uparrow}\_{ab} = \overline{\rho^{\uparrow\downarrow}\_{Ib,Ia}(0)} \\) です。

〔規約〕その他：

- 非スピン偏極の計算では \\( \rho \\) は 1 スピン分（[密度行列](density-mixing.md)）で、\\( n \\) も 1 スピン分です。エネルギーは 2 倍します（下記）。
- 軌道最適化（`orbitalOpt.Method`）を使った場合、U は縮約前の（原始）軌道に掛かります。占有行列は縮約係数で原始軌道の密度行列に戻してから作ります（TN6 の 6 節）。
- `full` は応力計算（`scf.stress.tensor`）と併用できず、入力の検査で停止します。非共線も応力と併用できません（[全エネルギーと力](energy-force.md)）。

### SCF の中での位置〔規約〕

- 占有行列とポテンシャル \\( v \\) は、電荷ミキシングの **後** の（入力側の）密度行列から作ります。この \\( v \\) が次の反復の Hamiltonian に入り、
  その反復で出力する `Uhub`・`Ucs`・`Uzs` も同じ \\( n \\) で評価します。
- 収束後に、最後の密度行列から占有行列を作り直してファイル `System.Name.DM_onsite` に書きます
  （原子・スピンごとの占有行列、その固有値と固有ベクトル。見出しは「Occupation Number in LDA+U and Constraint DFT」）。
- LDA+U・拘束・Zeeman 項のどれかがあると、再開時の最初の反復は保存した Hamiltonian を使います（[密度行列](density-mixing.md) の 6 節）。

## 2. 簡略形（`scf.DFTU.Type 1`）〔公開〕〔規約〕

〔公開〕TN6 Eq. 2–4、TN2 Eq. 101（Dudarev らの回転不変形）：

$$
E\_U = \frac12 \sum\_s U\_s \sum\_\sigma \left[ \mathrm{Tr}\\, n^{\sigma\sigma}\_s - \mathrm{Tr}\\, (n^{\sigma\sigma}\_s n^{\sigma\sigma}\_s) \right]
\quad\text{（共線）},\qquad
E\_U = \frac12 \sum\_s U\_s \left[ \mathrm{Tr}\\, n\_s - \mathrm{Tr}\\, (n\_s n\_s) \right]
\quad\text{（非共線、}2(2l+1)\text{ 次元）}
$$

$$
v^{\sigma\sigma'}\_{s,ab} = U\_s \left( \tfrac12 \delta\_{\sigma\sigma'}\delta\_{ab} - n^{\sigma'\sigma}\_{s,ba} \right)
$$

- 和はシェル内だけで、異なるシェル（異なる \\( l \\) や \\( \zeta \\)）の間の項はありません。入力の \\( U \\) は Dudarev 形の
  \\( U\_{\mathrm{eff}} = U - J \\) の意味で使われます（J は読みません）。
- 〔規約〕非スピン偏極では、1 スピン分の式を 2 倍したものを `Uhub` とします（ポテンシャルは 1 スピン分の式のまま）。
  〔公開〕マニュアル 40.1 節は DFT+U に `scf.SpinPolarization` の `on` か `nc` を求めていますが、既存実装は簡略形の非偏極計算を止めずに実行します。
- 冪等な（固有値が 0 か 1 の）\\( n\_s \\) では \\( E\_U = 0 \\) です。

## 3. 一般形（`scf.DFTU.Type 2`）〔公開〕〔規約〕

### 3.1 エネルギーとポテンシャル

〔公開〕Liechtenstein らの形と、Ryee–Han（Sci. Rep. 2018、JPCM 2018）の二重計数・非共線への拡張です。
シェルごとに Hartree–Fock 型の相互作用

$$
E\_{\mathrm{HF}}[n] = \frac12 \sum\_{abcd} U\_{abcd}\\, \mathcal N\_{ac}\\, \mathcal N\_{bd}
{}- \frac12 \sum\_{abcd} U\_{abdc} \sum\_{\sigma\sigma'} n^{\sigma\sigma'}\_{ac}\\, n^{\sigma'\sigma}\_{bd},\qquad
\mathcal N = n^{\uparrow\uparrow} + n^{\downarrow\downarrow}
$$

を作り（全体の実部をとる）、二重計数補正 \\( E\_{\mathrm{dc}} \\) を引きます。共線では \\( n^{\uparrow\downarrow} = 0 \\) とおけば

$$
E\_{\mathrm{HF}} = \frac12\sum U\_{abcd}\left(n^{\uparrow}\_{ac} n^{\downarrow}\_{bd} + n^{\downarrow}\_{ac} n^{\uparrow}\_{bd}\right)
{}+ \frac12\sum \left(U\_{abcd} - U\_{abdc}\right)\sum\_\sigma n^{\sigma}\_{ac} n^{\sigma}\_{bd}
$$

で、〔規約〕既存実装の共線の式と同じです。ポテンシャルは

$$
v^{\sigma\sigma'}\_{ac} = \delta\_{\sigma\sigma'} \sum\_{bd} U\_{abcd}\\, \mathcal N\_{bd} - \sum\_{bd} U\_{abdc}\\, n^{\sigma'\sigma}\_{bd} - v^{\mathrm{dc},\sigma\sigma'}\_{ac}
$$

です（\\( U\_{abcd} = U\_{badc} \\) を使って \\( E \\) を微分したもの）。

### 3.2 二重計数補正（`scf.dc.Type`）〔公開〕〔規約〕

| `scf.dc.Type` | \\( E\_{\mathrm{dc}} \\) | ポテンシャルの補正 |
|---|---|---|
| `sFLL`（既定） | \\( \tfrac12 U N(N-1) - \tfrac12 J \left[ \sum\_{\sigma\sigma'} N^{\sigma\sigma'} N^{\sigma'\sigma} - N \right] \\) | \\( v^{\mathrm{dc},\sigma\sigma'}\_{ac} = \delta\_{ac}\left[ U(N-\tfrac12)\delta\_{\sigma\sigma'} - J(N^{\sigma'\sigma} - \tfrac12\delta\_{\sigma\sigma'}) \right] \\) |
| `cFLL` | \\( \tfrac12 U N(N-1) - \tfrac14 J N(N-2) \\) | \\( v^{\mathrm{dc},\sigma\sigma'}\_{ac} = \delta\_{ac}\delta\_{\sigma\sigma'}\left[ U(N-\tfrac12) - J(\tfrac12 N - \tfrac12) \right] \\) |
| `sAMF` | なし（\\( E = E\_{\mathrm{HF}}[\tilde n] \\)） | なし（\\( v = v\_{\mathrm{HF}}[\tilde n] \\)） |
| `cAMF` | なし（\\( E = E\_{\mathrm{HF}}[\tilde n] \\)） | なし（\\( v = v\_{\mathrm{HF}}[\tilde n] \\)） |

ここで \\( N, N^{\sigma\sigma'} \\) はシェルの値です。AMF（around mean field）は占有行列から平均を引いた \\( \tilde n \\) を使います。

- `sAMF`：\\( \tilde n^{\sigma\sigma'} = n^{\sigma\sigma'} - \dfrac{N^{\sigma\sigma'}}{2l+1} I \\)（非共線ではスピン非対角成分にも同じ形で引く）。
- `cAMF`：\\( \tilde n^{\sigma\sigma} = n^{\sigma\sigma} - \dfrac{N}{2(2l+1)} I \\)、スピン非対角成分はそのまま。

共線では \\( \sum\_{\sigma\sigma'} N^{\sigma\sigma'}N^{\sigma'\sigma} = (N^\uparrow)^2 + (N^\downarrow)^2 \\) なので、sFLL は通常の
\\( \tfrac12 UN(N-1) - \tfrac12 J\sum\_\sigma N^\sigma(N^\sigma-1) \\) になります。
AMF で \\( v = v\_{\mathrm{HF}}[\tilde n] \\) が \\( E\_{\mathrm{HF}}[\tilde n(n)] \\) の正確な微分になるのは、回転不変なテンソルでは
\\( \mathrm{Tr}\\, v\_{\mathrm{HF}}[\tilde n] \\) の、引いた平均に対応する成分が 0 になるためです（検算で確認）。

- 〔検算〕\\( J = 0 \\) の sFLL は、\\( U = F^0 \\) の簡略形とエネルギー・ポテンシャルとも一致します。sFLL・sAMF は全体のスピン回転で不変です。
- 〔公開〕マニュアル 40.1.2：先頭の `s`・`c` は、LDA/GGA 部分をスピン密度汎関数（s）で扱うか電荷密度汎関数（c）で扱うかの区別です。
  〔規約〕`cFLL`・`cAMF` では、SCF の **2 反復目以降**、格子上の交換相関を両スピンの密度を平均した密度（\\( n\_\uparrow = n\_\downarrow = n/2 \\)）で評価します
  （1 反復目はスピン分極した XC）。簡略形（Type 1）は暗黙に sFLL 相当です。
- 〔公開〕軌道分極の促進（4 節）は `cFLL`・`cAMF` との併用を推奨しない、とマニュアルにあります。

### 3.3 Coulomb テンソル〔公開〕〔規約〕

〔公開〕複素球面調和関数 \\( Y\_{lm} \\)（Condon–Shortley）での行列要素は、Slater 積分 \\( F^k \\) と Gaunt 係数で

$$
\langle m\_1 m\_2 \vert V \vert m\_3 m\_4\rangle = \delta\_{m\_1+m\_2,\\, m\_3+m\_4} \sum\_{k=0,2,\dots,2l} c^k(m\_1, m\_3)\\, c^k(m\_4, m\_2)\\, F^k,
$$

$$
c^k(m, m') = \sqrt{\tfrac{4\pi}{2k+1}} \int d\Omega\\, Y^{\ast}\_{lm} Y\_{k,m-m'} Y\_{lm'}
{}= (-1)^{m} (2l+1) \begin{pmatrix} l & k & l\\ 0&0&0\end{pmatrix}\begin{pmatrix} l & k & l\\ -m & m-m' & m'\end{pmatrix}
$$

です。実球面調和関数 \\( S\_a = \sum\_m T\_{am} Y\_{lm} \\)（[PAO 基底](pao-basis.md) の \\( U^{(l)} \\)）へは
\\( U\_{abcd} = \sum \overline{T\_{a m\_1}}\\,\overline{T\_{b m\_2}}\\, T\_{c m\_3} T\_{d m\_4} \langle m\_1 m\_2 \vert V \vert m\_3 m\_4\rangle \\) で変換し、実数になります。
平均値は \\( U = (2l+1)^{-2}\sum\_{ab} U\_{abab} = F^0 \\)、\\( U - J = [2l(2l+1)]^{-1}\sum\_{ab}(U\_{abab} - U\_{abba}) \\) です。

〔規約〕\\( F^k \\) は入力の \\( U, J \\) から次のように決めます（`scf.Yukawa off`、既定）。

| \\( l \\) | \\( F^0 \\) | \\( F^2 \\) | \\( F^4 \\) | \\( F^6 \\) |
|---|---|---|---|---|
| 0 | \\( U \\) | — | — | — |
| 1 | \\( U \\) | \\( 5J \\) | — | — |
| 2 | \\( U \\) | \\( 14J/(1+r) \\) | \\( rF^2 \\) | — |
| 3 | \\( U \\) | \\( 6435J/(286 + 195\cdot 0.668 + 250\cdot 0.494) \\) | \\( 0.668\\, F^2 \\) | \\( 0.494\\, F^2 \\) |

\\( r \\) は `scf.Slater.Ratio`（既定 0.625）です。f 軌道の比 0.668・0.494 は固定で、入力では変えられません。
いずれも \\( J = (F^2+F^4)/14 \\)（d）、\\( F^2/5 \\)（p）、\\( (286F^2 + 195F^4 + 250F^6)/6435 \\)（f）を満たします（検算）。

〔規約〕既存実装の実球面調和関数への変換は、[PAO 基底](pao-basis.md) の並びと符号（すべて正の係数）と一致しています。
\\( U\_{abcd} \\) は個々の実軌道の符号に依存するので、再実装でも PAO 基底と同じ \\( T \\) を使います。

### 3.4 Yukawa 型の遮蔽による J の推定（`scf.Yukawa on`）〔公開〕〔規約〕

〔公開〕マニュアル 40.2.2：入力の \\( U \\) だけを使い、遮蔽 Coulomb ポテンシャル \\( e^{-\lambda r}/r \\) の遮蔽長 \\( \lambda \\) を
\\( F^0 = U \\) となるように決め、同じ \\( \lambda \\) で \\( F^2, F^4, F^6 \\) と \\( J \\) を作ります。`Hund.J.values` と `scf.Slater.Ratio` は無視されます。

$$
F^k(\lambda) = (2k+1)\lambda \int\_0^\infty\\!\\!\int\_0^\infty dr\\, dr'\\, r^2 r'^2 R^2(r) R^2(r')\\, i\_k(\lambda r\_<)\\, \tilde k\_k(\lambda r\_>)
$$

\\( R \\) はそのシェルの PAO 動径関数、\\( i\_k \\) は変形球 Bessel 関数（\\( i\_0(x) = \sinh x/x \\)）、\\( \tilde k\_k \\) は \\( \tilde k\_0(x) = e^{-x}/x \\) と同じ漸化式を満たす第 2 種の関数です
（\\( e^{-\lambda R}/R = \lambda\sum\_k (2k+1) i\_k(\lambda r\_<)\tilde k\_k(\lambda r\_>) P\_k(\cos\gamma) \\)、\\( \lambda\to 0 \\) で通常の \\( r\_<^k/r\_>^{k+1} \\)。検算）。

〔規約〕数値の決め方（結果に影響するもの）：

- \\( \lambda \\) は区間 \\( [0.1, \lambda\_{\max}] \\) bohr⁻¹ の二分法で求めます。\\( \lambda\_{\max} \\) は \\( U \ge 1.5 \\) Hartree なら 10、\\( U \ge 0.5 \\) Hartree なら 20、それ以外は 40 です。
  \\( \lvert F^0 - U\rvert \le 0.0005 \\) eV で止め、200 回で収束しなければ警告を出してそのまま進みます。
- 動径積分は \\( r \in [0, r\_c + 1] \\) bohr（\\( r\_c \\) は PAO のカットオフ半径）の Simpson 則で、外側の刻み 0.01 bohr、内側の累積積分の刻み 0.001 bohr です。
- \\( J \\) は上の表の関係（p：\\( F^2/5 \\)、d：\\( (F^2+F^4)/14 \\)、f：\\( (286F^2+195F^4+250F^6)/6435 \\)）から求め、出力に \\( \lambda \\)、\\( U \\)、\\( J \\)、\\( F^k \\)、\\( F^4/F^2 \\) を書きます。

### 3.5 一般形の適用範囲〔規約〕

- テンソルを作るのは \\( U > 10^{-5} \\) Hartree または \\( J > 10^{-5} \\) Hartree のシェルだけで、それ以外のシェルには何も加えません。
- **既存実装の一般形は \\( \zeta = 0, 1 \\) のシェルだけを正しく扱います。** \\( \zeta \ge 2 \\) のシェルに 0 でない U・J を与えると、シェルの先頭位置の計算が正しくなく、結果は未定義です。
- 非スピン偏極の計算での一般形は、既存実装では定義されていません（マニュアルもスピン偏極を求めています）。
- 一般形は `scf.Hubbard.U on` のときだけ有効で、`scf.DFTU.Type` の既定は 1 です。

## 4. 軌道分極の促進〔公開〕〔規約〕

〔公開〕マニュアル 40.1.3、TN6 の 5 節：SCF の最初の数反復で、占有行列を軌道分極した形に置き換え、
軌道分極した解に入りやすくします。入力は `Atoms.SpeciesAndCoordinates` の最後の欄のフラグです（`on`・`off`）。

〔規約〕既存実装の規則：

- 条件：`scf.Hubbard.U on`、SCF 反復の番号が 1〜8、再開ファイルから密度行列を読み込んでいない（つまり MD・構造最適化の 2 ステップ目以降は働かない）、
  その原子のフラグが `off` でない。
- 対象：その原子の \\( l \ge 2 \\) の各シェルのうち **\\( \zeta = 0 \\) のものだけ**。行列は、`scf.Hubbard.Occupation` の選択によらず
  onsite の行列 \\( \rho\_{II}(0) \\) の該当部分です。
- 手順：行列（下の各場合）を対角化し、固有値の和 \\( N \\) が 0.3 を超えるときだけ、固有値の大きい順に \\( \lfloor N\rfloor \\) 個の固有ベクトルに占有 1、
  次の 1 個に端数 \\( N - \lfloor N\rfloor \\)、残りに 0 を与えて行列を作り直し、その反復の占有行列（ポテンシャルとエネルギーに使うもの）の該当ブロックを置き換えます。
  電子数は保存されます。
  - 非スピン偏極：\\( (2l+1) \\) 次元の 1 スピン分。フラグ `ex` では、端数を与えた準位と、占有 1 の準位のうち最も低いものの占有を入れ替えます（第一励起配置）。
  - 共線スピン偏極：両スピンのブロックを並べた \\( 2(2l+1) \\) 次元の行列（スピン非対角は 0）を対角化し、スピン対角ブロックだけを戻します。
    **端数は「次の準位」ではなく、さらに 1 つ下の準位に与えられます**（既存実装の挙動。占有 0 の準位が 1 つ挟まる）。
  - 非共線：\\( 2(2l+1) \\) 次元の複素行列。端数は次の準位です。作り直した行列 \\( P \\) を、各スピンブロックの軌道添字について転置して
    （\\( n^{\sigma\sigma'}\_{ba} = P^{\sigma\sigma'}\_{ab} \\)）格納します。
  - 共線スピン偏極と非共線でのフラグ `ex` は未対応で、該当ブロックの占有がすべて 0 になります。
- 既存実装は対角化の前に、行列要素へ \\( 10^{-13} \\)〜\\( 10^{-12} \\) 程度の乱数の摂動を加えます。結果への影響は無視できるので、再実装では不要です〔設計〕。

## 5. 非共線スピンの向きの拘束（`scf.Constraint.NC.Spin`）〔公開〕〔規約〕

〔公開〕マニュアル 41 節：各原子のスピンの向きを初期の向きに（`on2` では大きさも初期値に）近づけるペナルティ汎関数を加えます。
強さは `scf.Constraint.NC.Spin.V`（eV、既定 0）で、拘束する原子は座標の欄の 12 列目を 1 にします。

〔規約〕原子 \\( I \\) の 2×2 スピン密度行列を、その原子の全軌道の dual 占有行列の対角和から

$$
T^{(I)}\_{\sigma\sigma'} = \sum\_{\alpha \in I} n^{\sigma'\sigma}\_{I,\alpha\alpha}
$$

と作ります（電荷密度の規約、TN2 Eq. 3 の \\( n\_{\sigma\sigma'} \\) と同じ向きで、\\( n \\) の転置にあたる）。全電子数と磁化は

$$
N\_I = \mathrm{Tr}\\, T,\qquad m\_{I,c} = \mathrm{Tr}\\,(\sigma\_c T) \quad\Leftrightarrow\quad
\mathbf m = (2\\,\mathrm{Re}\\, T\_{\uparrow\downarrow},\ -2\\,\mathrm{Im}\\, T\_{\uparrow\downarrow},\ T\_{\uparrow\uparrow} - T\_{\downarrow\downarrow})
$$

で、[グリッド](grid-poisson-xc.md) の 2.4 節と同じ規約です。目標の行列は、初期方向 \\( \hat{\mathbf e}\_0 = \hat{\mathbf e}(\theta\_0, \varphi\_0) \\)（下の 8 節）を使って

$$
T\_0 = \tfrac12\left( N\_I\\, \sigma\_0 + M\\, \hat{\mathbf e}\_0\cdot\boldsymbol\sigma \right),\qquad
M = \begin{cases} \lVert \mathbf m\_I \rVert & (\texttt{on}) \\ M\_0 = \lvert N^{\uparrow}\_{I,0} - N^{\downarrow}\_{I,0}\rvert & (\texttt{on2}) \end{cases}
$$

です（\\( N^{\uparrow,\downarrow}\_{I,0} \\) は入力の初期占有）。TN2 Eq. 16 の回転行列 \\( U(\theta\_0,\varphi\_0) \\) で
\\( T\_0 = U^\dagger \mathrm{diag}(\lambda\_+, \lambda\_-) U \\)（\\( \lambda\_\pm = (N\_I \pm M)/2 \\)）とも書けます（検算）。ペナルティは Frobenius ノルムで

$$
E\_{\mathrm{cs}} = V \sum\_{I:\\,\text{拘束}} \lVert T^{(I)} - T\_0^{(I)} \rVert\_F^2
{}= \frac{V}{2}\sum\_I \lVert \mathbf m\_I - M\hat{\mathbf e}\_0 \rVert^2
$$

で、`on` では \\( E\_{\mathrm{cs}} = V\lVert\mathbf m\rVert^2 (1 - \cos\gamma) \\)（\\( \gamma \\) は \\( \mathbf m \\) と \\( \hat{\mathbf e}\_0 \\) のなす角）、`on2` では \\( \tfrac V2 \lVert \mathbf m - M\_0\hat{\mathbf e}\_0\rVert^2 \\) です。
出力の `Ucs` はこの値です。ポテンシャルは、有効磁場 \\( \mathbf h = \partial E\_{\mathrm{cs}}/\partial \mathbf m \\) を使って

$$
v^{\sigma\sigma'}\_{I,\alpha\beta} \mathrel{+}= \delta\_{\alpha\beta}\\, (\mathbf h\cdot\boldsymbol\sigma)\_{\sigma\sigma'}\quad(\alpha \in I),\qquad
\mathbf h = \begin{cases} V\left( 2\mathbf m - \hat{\mathbf m}\\,(\mathbf m\cdot\hat{\mathbf e}\_0) - \lVert\mathbf m\rVert\\,\hat{\mathbf e}\_0 \right) & (\texttt{on}) \\ V(\mathbf m - M\_0\hat{\mathbf e}\_0) & (\texttt{on2}) \end{cases}
$$

で、これは \\( E\_{\mathrm{cs}} \\) の \\( n \\) についての正確な微分です（検算）。\\( v \\) は dual の Hamiltonian（6 節）として加わります。

- 〔規約〕既存実装は、目標 \\( T\_0 \\) の \\( T \\) への依存（`on` では固有値を通して、`on2` では \\( N\_I \\) を通して）を連鎖律で含めたポテンシャルを使います。
  その微分が上の正確な勾配と丸め誤差の範囲で一致するかは確認していません（〔未確定〕）。
- 〔公開〕マニュアル 41 節：拘束は物理的な外場ではなく数値的な制約で、DFT+U と併用できます。局所解を避けるため、容易軸から 30° 程度ずつ回し、前の結果から再開する使い方が推奨されています。

## 6. Zeeman 項（`scf.NC.Zeeman.Spin`・`scf.NC.Zeeman.Orbital`）〔公開〕〔規約〕

〔公開〕TN2 Eq. 139–150、マニュアル 43 節：一様磁場とスピン・軌道磁気モーメントの相互作用 \\( E\_z = \sum\_I \mathbf B^s\_I\cdot\mathbf s\_I + \tfrac12 \sum\_I \mathbf B^o\_I\cdot\mathbf l\_I \\)。
磁場の大きさは `scf.NC.Mag.Field.Spin`・`scf.NC.Mag.Field.Orbital`（Tesla、既定 0）、向きは座標の欄の 8・9 列（スピン）、10・11 列（軌道）の角度です。

〔規約〕

- 単位：\\( B[\text{a.u.}] = B[\mathrm T] / 2.35051742\times 10^{5} \\)（CODATA 2018 の \\( 2.35051757\times10^5 \\) T とは相対 \\( 6\times10^{-8} \\) の差）。
- **12 列目が 1 の原子にだけ** 働きます（スピン・軌道とも）。向きはその原子の初期方向（8 節の補正後）です。
- スピン：\\( \mathbf s\_I = \tfrac12 \mathbf m\_I \\)（5 節の \\( T^{(I)} \\) から）、\\( \mathbf b\_I = -B\_s\\, \hat{\mathbf e}\_{s,I} \\) として

$$
U\_{\mathrm{zs}} = \sum\_I \mathbf b\_I\cdot \mathbf s\_I = -\tfrac12 B\_s \sum\_I \hat{\mathbf e}\_{s,I}\cdot\mathbf m\_I,\qquad
v^{\sigma\sigma'}\_{I,\alpha\beta} \mathrel{+}= \delta\_{\alpha\beta}\\, \tfrac12 (\mathbf b\_I\cdot\boldsymbol\sigma)\_{\sigma\sigma'}
$$

  正の \\( B\_s \\) は、OpenMX の表示するスピンモーメント（\\( \mathbf m \\)）を \\( \hat{\mathbf e}\_s \\) にそろえる向きに働きます。
  〔公開〕マニュアル 37 節のとおり、OpenMX の表示するモーメントは物理的な電子の磁気モーメントと逆符号の規約です。
- 軌道：\\( \mathbf a\_I = -\tfrac12 B\_o\\, \hat{\mathbf e}\_{o,I} \\)、\\( \mathbf L\_I \\) を 7 節の軌道モーメントとして

$$
U\_{\mathrm{zo}} = \sum\_I \mathbf a\_I\cdot\mathbf L\_I,\qquad
v^{\sigma\sigma}\_{I,\alpha\beta} \mathrel{+}= \mathbf a\_I\cdot\langle \phi\_{I\alpha}\vert \hat{\mathbf L}\vert \phi\_{I\beta}\rangle \quad(\text{両スピン、スピン対角のみ})
$$

  （\\( \langle\phi\_{I\alpha}\vert\hat L\_c\vert\phi\_{I\beta}\rangle \\) は純虚数、7 節）。この \\( v \\) も dual の Hamiltonian として加わるので、
  \\( \sum H\rho \\) への寄与は dual の \\( \mathbf a\cdot\mathbf L \\) と一致します。ただし \\( U\_{\mathrm{zo}} \\) の \\( \mathbf L\_I \\) は直前の対角化の出力密度行列から評価した値です（ミキシング前）。
- 入力の検査：いずれも非共線（`scf.SpinPolarization NC`）と `scf.Hubbard.Occupation dual` が必要です（拘束も同じ）。
  軌道の Zeeman 項はスピン軌道相互作用（`scf.SpinOrbit.Coupling on`）が必要です。スピンの Zeeman 項と `scf.Constraint.NC.Spin` は同時に使えません。
  条件を満たさないと停止します。
- `Uzs`・`Uzo` は `Utot` に含まれます（[全エネルギーと力](energy-force.md)）。外場なしのエネルギーと比べるときは差し引きます。

## 7. 軌道磁気モーメント〔公開〕〔規約〕

〔公開〕TN2 Eq. 151–162、マニュアル 39 節：非共線計算で、各原子の軌道角運動量の期待値を局在基底で評価します。スピン軌道相互作用がなければ 0 です。

〔規約〕原子 \\( I \\) 上の軌道角運動量の行列要素を

$$
\langle \phi\_{I\alpha} \vert \hat L\_c \vert \phi\_{I\beta}\rangle = i\\, \Lambda^{c}\_{\alpha\beta},\qquad
\Lambda^c\_{\alpha\beta} = \Big(\int r^2 R\_{\alpha} R\_{\beta}\\, dr\Big)\\, \mathrm{Im}\\,\langle S\_{l a}\vert \hat L\_c \vert S\_{l b}\rangle\\,\delta\_{l\_\alpha l\_\beta}
$$

（\\( \hbar = 1 \\)、実の \\( \Lambda \\) は反対称）とし、スピン対角成分の密度行列の虚部の dual 射影

$$
\tilde n^{\sigma}\_I = \tfrac12\left( \mathrm{Im}\\,\rho^{\sigma\sigma}\\, S + S\\, \mathrm{Im}\\,\rho^{\sigma\sigma}\right)\_{II}
$$

から

$$
L\_{I,c} = -\sum\_{\sigma}\sum\_{\alpha\beta\in I} \Lambda^c\_{\alpha\beta}\\, \tilde n^{\sigma}\_{I,\alpha\beta}
{}= \mathrm{Re}\sum\_\sigma\sum\_{\alpha\beta}\langle\phi\_{I\alpha}\vert\hat L\_c\vert\phi\_{I\beta}\rangle\\, n^{\sigma\sigma}\_{I,\alpha\beta}
$$

と定めます（\\( n \\) は dual。\\( \mathrm{Re}\\,n \\) は対称なので反対称の \\( \Lambda \\) とは打ち消し合う）。
例えば p 軌道（\\( x, y, z \\) の順）で \\( \langle x\vert \hat L\_z\vert y\rangle = -i \\)、\\( \langle y\vert\hat L\_z\vert x\rangle = i \\) です（TN2 Eq. 162 と同じ）。

- 単位は μB として、\\( \mathbf L\_I \\) をそのまま表示します（物理的な軌道磁気モーメントは \\( -\mu\_B \mathbf L \\)。表示は逆符号の規約）。
- 各 SCF 反復の対角化の後に、出力の（ミキシング前の）密度行列から計算します。全体の和 \\( \sum\_I \mathbf L\_I \\) とその角度を標準出力に、
  原子ごとの大きさと角度、軌道・多重度ごとの分解をファイル `System.Name.OM`（`System.Name.out` の「Orbital moments」の節にも含まれる）に書きます。
- 原子ごとの大きさと角度 \\( (\theta, \varphi) \\) は \\( \mathbf L\_I \\) の極座標です。

## 8. 初期磁気配置の入力〔公開〕〔規約〕

〔公開〕マニュアル 37 節：`Atoms.SpeciesAndCoordinates` の各行の欄は次のとおりです。

| 列 | 共線（`on`・`off`） | 非共線（`NC`） |
|---|---|---|
| 1 | 通し番号 | 通し番号 |
| 2 | 種 | 種 |
| 3–5 | 座標 | 座標 |
| 6 | 上向きスピンの初期占有 \\( N^\uparrow\_0 \\) | 同左 |
| 7 | 下向きスピンの初期占有 \\( N^\downarrow\_0 \\) | 同左 |
| 8 | 軌道分極フラグ | スピンの角度 \\( \theta\_s \\)（度） |
| 9 | — | スピンの角度 \\( \varphi\_s \\)（度） |
| 10 | — | 軌道（モーメント・磁場）の角度 \\( \theta\_o \\)（度） |
| 11 | — | 軌道の角度 \\( \varphi\_o \\)（度） |
| 12 | — | 拘束・Zeeman 項を働かせるフラグ（整数、1 で有効） |
| 13 | — | 軌道分極フラグ |

〔規約〕

- 初期密度：原子密度を \\( N^\uparrow\_0 + N^\downarrow\_0 \\) で拡大縮小し、上向き・下向きに割り振ります（[密度行列](density-mixing.md) の 6 節）。
  非共線では方向 \\( \hat{\mathbf e}(\theta\_s, \varphi\_s) \\) に磁化を向けた 2×2 の密度を作ります。
- **多数スピンの向きへの付け替え**：非共線で \\( N^\uparrow\_0 < N^\downarrow\_0 \\) の行は、\\( N^\uparrow\_0 \\) と \\( N^\downarrow\_0 \\) を入れ替え、
  スピンと軌道の方向をそれぞれ反対向き（\\( \hat{\mathbf e} \to -\hat{\mathbf e} \\)）に付け替えます。したがって内部の「初期方向」は常に多数スピンの磁化の向きで、
  拘束（5 節）・Zeeman 項（6 節）の \\( \hat{\mathbf e}\_0 \\) もこの付け替え後の方向です。
- 角度の小さな揺らぎ：\\( \lvert\theta\rvert < 10^{-10} \\) 度のとき、\\( \theta\_s \\)・\\( \theta\_o \\) に \\( 10^{-5} \\) 度未満の乱数を加えます。
  さらに \\( \theta\_s, \varphi\_s \\) には \\( 10^{-10} \\) rad 未満の乱数を加えます（方位角の特異点を避けるため）。結果への影響は無視できます。
- 軌道分極フラグは `off`・`on`・`ex`（大文字・小文字を区別しない）で、`scf.Hubbard.U on` のときだけ読みます。ほかの文字列は停止します。
  〔公開〕マニュアルは省略時を `off` としていますが、既存実装は欄がない行で値を初期化しません（直前の行の値が残る）。全原子に明示するのが安全です。
- 欄が足りない行では、既存実装は足りない欄の値を更新しません（角度・フラグの配列の初期値が残る）。〔設計〕再実装では角度 0、フラグ 0、軌道分極 `off` を既定値とします。

## 9. DFT+U の入力〔公開〕〔規約〕

| キーワード | 値 | 既定 | 意味 |
|---|---|---|---|
| `scf.Hubbard.U` | on／off | off | DFT+U を使う |
| `scf.DFTU.Type` | 1／2 | 1 | 簡略形／一般形 |
| `scf.Hubbard.Occupation` | onsite／full／dual | dual | 占有行列の定義（拘束・Zeeman 項でも使う） |
| `<Hubbard.U.values` | ブロック | — | シェルごとの U（eV） |
| `<Hund.J.values` | ブロック | — | シェルごとの J（eV、Type 2 で `scf.Yukawa off` のとき） |
| `scf.dc.Type` | sFLL／sAMF／cFLL／cAMF | sFLL | 二重計数補正（Type 2） |
| `scf.Slater.Ratio` | 実数 | 0.625 | d の \\( F^4/F^2 \\)（Type 2） |
| `scf.Yukawa` | on／off | off | 遮蔽長から J と \\( F^k \\) を推定（Type 2） |
| `scf.Constraint.NC.Spin` | off／on／on2 | off | スピンの向き（と大きさ）の拘束 |
| `scf.Constraint.NC.Spin.V` | 実数（eV） | 0 | 拘束の強さ |
| `scf.NC.Zeeman.Spin`、`scf.NC.Zeeman.Orbital` | on／off | off | Zeeman 項 |
| `scf.NC.Mag.Field.Spin`、`scf.NC.Mag.Field.Orbital` | 実数（T） | 0 | 磁場の大きさ |

〔規約〕

- 選択肢の文字列は大文字・小文字を区別しません。
- `<Hubbard.U.values` の各行は、種の名前のあとに「ラベル 値」の組を、その種の基底の \\( (l, \zeta) \\) について
  \\( l = 0, 1, \dots \\)、各 \\( l \\) の中で \\( \zeta = 0, 1, \dots \\) の順に、**基底の指定（例：`s2p2d2f1`）と同じ数だけ** 並べます
  （例：`Fe 1s 0.0 2s 0.0 1p 0.0 2p 0.0 1d 5.0 2d 0.0 1f 0.0`）。ラベル（`1d` など）は読み飛ばされ、照合されません。
  順番と数が合っていないと値が別のシェルに割り当てられます。行は種の数だけ読み、種の順序は任意です。`<Hund.J.values` も同じ形式です。
- eV から Hartree への換算は 27.2113845 で割ります。
- `scf.DFTU.Type` 1 では、`scf.dc.Type` と J は読みません。

## 10. Hamiltonian への組み込み〔公開〕〔規約〕

〔公開〕TN6 Eq. 21–28、TN2 Eq. 122–124：シェルのポテンシャル \\( v\_s \\) を原子ごとのブロック対角行列 \\( V \\)（\\( V\_{II} \\) の中のシェル \\( s \\) の部分が \\( v\_s \\)、
シェル間と原子間は 0）にまとめ、射影の種類に応じて

| 占有行列 | Hamiltonian への寄与 \\( H^U \\) | 原子対の範囲 |
|---|---|---|
| onsite | \\( H^U\_{II}(0) = V\_{II} \\) | 原子自身だけ |
| full | \\( H^U = S V S \\) | \\( S \\) の 2 段の到達範囲（既存の原子対の集合に含まれる成分だけ） |
| dual | \\( H^U = \tfrac12 (V S + S V) \\)、すなわち \\( H^U\_{I\alpha,J\beta}(\mathbf L) = \tfrac12[ (V\_{II} S\_{IJ}(\mathbf L))\_{\alpha\beta} + (S\_{IJ}(\mathbf L) V\_{JJ})\_{\alpha\beta} ] \\) | \\( S \\) と同じ |

を加えます。非共線では \\( V \\) の各スピンブロック \\( V^{\sigma\sigma'} \\) に同じ式を当て、虚部は Hamiltonian の虚部（スピン軌道相互作用の項と同じ場所）に足します。
拘束・Zeeman 項（5・6 節）の \\( v \\) は同じ \\( V \\) に足され、`scf.Hubbard.Occupation` の射影で Hamiltonian に入ります（入力の検査で dual に限られる）。

- 定義から \\( \sum\_{ab} H^U\_{ab}\rho\_{ab} = \sum\_s\sum\_{ab} v\_{s,ab}\\, n\_{s,ab} \\) が成り立ちます（検算。共線の非偏極では両辺とも 1 スピン分）。
- 〔規約〕full の \\( S V S \\) は、\\( (I, J, \mathbf L) \\) が疎ブロック形式の原子対の集合に含まれる成分だけを作ります（範囲外は捨てる）。

## 11. 力と応力〔公開〕〔規約〕

〔公開〕TN6 Eq. 29–30：\\( \partial E\_U/\partial\tau = \sum v\\, \partial n/\partial \tau \\)。\\( n \\) の \\( \rho \\) を通した変化は、
Hamiltonian の \\( H^U \\) とエネルギー重み付き密度行列の \\( \partial S \\) の項で既に含まれるので（[全エネルギーと力](energy-force.md) の Pulay の恒等式）、
残るのは \\( \rho \\) を固定したときの \\( S \\) を通した変化です。\\( V \\) を固定して

$$
\frac{\partial E\_U}{\partial \mathbf R\_K}\bigg|\_{\rho} = g\_s \sum\_{\sigma\sigma'} \sum\_{\mathrm{pairs}} \rho^{\sigma\sigma'}\_{ab}\\, \frac{\partial H^{U,\sigma\sigma'}\_{ab}}{\partial \mathbf R\_K}\bigg|\_{V},\qquad
\frac{\partial H^U}{\partial \mathbf R\_K} = \begin{cases} 0 & \text{onsite}\\ \partial S\\, V S + S V \partial S & \text{full}\\ \tfrac12(V\\,\partial S + \partial S\\, V) & \text{dual}\end{cases}
$$

です（\\( g\_s \\) は非偏極で 2、それ以外 1。実部をとる）。\\( \partial S/\partial\mathbf R\_K \\) は原子 \\( K \\) を含む原子対ブロックだけが 0 でなく、
重なり行列の微分（[二中心積分](two-center-overlap-kinetic.md)）をそのまま使います。力は \\( \mathbf F\_K = -\partial E/\partial\mathbf R\_K \\) です。
有限差分で確認しています（`tools/check_dftu_magnetism.py`、dual・full、共線・非共線）。

〔規約〕

- 拘束・Zeeman 項も同じ \\( V \\) に入るので、同じ式で力に含まれます。軌道の Zeeman 項の \\( \Lambda \\) は原子自身のブロックだけで位置に依存しません。
- 既存実装は、full と dual のそれぞれについて、共線・非共線とも上の項を加えます（onsite は 0）。
- 応力：full は応力と併用できません（入力で停止）。非共線は応力計算自体が未対応です。共線の dual の応力の寄与は〔未確定〕です。

## 12. 磁気異方性エネルギー（第二変分法）〔公開〕〔規約〕

〔公開〕マニュアル 42 節：強磁性体について、(1) スピン軌道相互作用なしの共線計算を収束させ、(2) その再開ファイルを
`scf.restart c2n` で非共線計算に変換し、スピンの向き `scf.Restart.Spin.Angle.Theta`・`scf.Restart.Spin.Angle.Phi`（度）でスピン軌道相互作用を含めて
一度だけ対角化し、(3) 向きごとの `Uele`（バンドエネルギー）の差を MAE とします（力の定理。入力密度だけに依存する二重計数項は向きによらないので `Utot` ではなく `Uele` を比べる）。

〔規約〕

- `scf.Restart.Spin.Angle.Theta`・`...Phi` の既定は 0 です。
- 変換では、格子上の共線の密度 \\( n\_\uparrow, n\_\downarrow \\) を、全格子点で同じ向き \\( \hat{\mathbf e}(\theta,\varphi) \\) に回して
  \\( \underline n = U^\dagger \mathrm{diag}(n\_\uparrow, n\_\downarrow) U = \tfrac12[(n\_\uparrow+n\_\downarrow)\sigma\_0 + (n\_\uparrow - n\_\downarrow)\hat{\mathbf e}\cdot\boldsymbol\sigma] \\)
  （\\( U \\) は TN2 Eq. 16 の行列、電荷密度の規約）とします。
- 変換と再開の詳細（密度行列・Hamiltonian の扱い、`scf.maxIter` との関係）は〔未確定〕です。

## エッジケース〔規約〕〔設計〕

- U がすべて 0 のシェルは、簡略形では何も加えず（\\( v = 0 \\)）、一般形では閾値 \\( 10^{-5} \\) Hartree 未満なら扱いません。
- 一般形の \\( \zeta \ge 2 \\) のシェル、非偏極の一般形：既存実装では未定義（3.5 節）。〔設計〕再実装は入力の段階でエラーにします。
- 軌道分極の促進の共線スピン偏極での端数の割り当て、`ex` の未対応（4 節）は既存実装の挙動です。〔設計〕照合のときだけ再現できるようにし、既定は「端数は次の準位、`ex` は非偏極以外ではエラー」とします。
- 拘束の `on` で \\( \lVert\mathbf m\rVert \to 0 \\) の原子は \\( \hat{\mathbf m} \\) が定まりません。\\( E\_{\mathrm{cs}} \\) と \\( \mathbf h \\) は \\( \lVert\mathbf m\rVert \\) について連続に 0 になるので、\\( \lVert\mathbf m\rVert < 10^{-12} \\) では \\( \mathbf h = 0 \\) とします〔設計〕。
- 占有行列はミキシング後の密度行列から作るので、SCF の途中では固有値が \\( [0,1] \\) を外れることがあります。エネルギー・ポテンシャルの式はそのまま使います（切り詰めない）。
- onsite・full では \\( \sum \mathrm{Tr}\\,n \neq N\_{\mathrm e} \\) で、シェルの占有が 1 を超えることがあります。

## 実装への要求〔設計〕

- **占有行列の射影**：`Projector`（onsite・full・dual）を列挙型にし、(a) \\( \rho \\) と \\( S \\) から各シェル（または各原子）の \\( n \\) を作る関数、
  (b) \\( V \\) と \\( S \\) から \\( H^U \\) を作る関数、(c) \\( V \\)・\\( \rho \\)・\\( \partial S \\) から力の寄与を作る関数、の 3 つを持たせます。
  3 つは同じ規約（\\( \sum H^U\rho = \sum v n \\)）を満たすことをテストで保証します。
- **相互作用汎関数**：シェルの \\( n\_s \\) を受け取り \\( (E\_s, v\_s) \\) を返すインターフェースにします。実装は `Dudarev { u }`、
  `General { tensor, u, j, dc: {SFLL, SAMF, CFLL, CAMF} }`。共線は非共線の特別な場合（スピン非対角 0）として同じコードで扱い、非偏極は呼び出し側で \\( g\_s = 2 \\) を掛けます。
- **Coulomb テンソル**：\\( (l, F^k) \\) から \\( U\_{abcd} \\) を作る純関数と、\\( (U, J, r) \\) から \\( F^k \\) を作る関数、Yukawa 推定を分けます。実球面調和関数の変換行列は PAO 基底のものを共有します。
- **原子ごとのスピン項**：拘束（`on`・`on2`）と Zeeman 項を、原子の \\( T^{(I)} \\)（と \\( \mathbf L\_I \\)）を受け取り \\( (E, \mathbf h) \\) を返す項として実装し、
  \\( v \\) への変換（\\( \mathbf h\cdot\boldsymbol\sigma \\) を原子の全軌道の対角に置く）は共通にします。拘束は正確な勾配を使います。
- **SCF との接続**：ミキシング後の \\( \rho \\) → \\( n \\) → \\( (E, V) \\) → 次の反復の \\( H \\)、の順を明示し、エネルギーと Hamiltonian が同じ \\( n \\) から作られることを型で保証します。
  `Uhub`・`Ucs`・`Uzs`・`Uzo` を別々に返します。
- **入力の検査**：8・9 節の条件（非共線・dual・スピン軌道相互作用の要否、拘束とスピン Zeeman の排他、full と応力、`Hubbard.U.values` の個数と順序）を
  読み込み時にすべて検査し、ラベル（`1d` など）も基底と照合します（既存実装より厳しくする）。
- **出力**：`System.Name.DM_onsite` 相当（原子・スピンごとの占有行列とその固有値・固有ベクトル）と、軌道モーメントの分解を出力できるようにします。
- **既存実装との照合用の切り替え**：軌道分極の促進の端数の位置と非共線での転置、`ex` の扱いを、照合モードでだけ再現できるようにします。

## テスト

### 解析的・不変量（`tools/check_dftu_magnetism.py` で確認済みのもの）

| 内容 | 許容誤差 |
|---|---|
| 実球面調和関数（p, d）が PAO 基底の表と一致 | \\( 10^{-14} \\) |
| \\( U\_{abcd} \\) が実数、対称性 \\( U\_{abcd} = U\_{badc} = U\_{cdab} = U\_{cbad} \\) | \\( 10^{-14} \\) |
| \\( (2l+1)^{-2}\sum U\_{abab} = F^0 \\)、平均から求めた \\( J \\) が入力の \\( J \\)（l = 1, 2, 3） | \\( 10^{-14} \\) |
| 任意の回転 \\( D \\) で \\( U\_{abcd} \\) が不変 | \\( 10^{-11} \\) |
| 簡略形・一般形（4 種の二重計数）、共線・非共線で \\( v = \partial E/\partial n \\)（中心差分、相対） | \\( 10^{-8} \\) |
| \\( J = 0 \\) の sFLL ＝ \\( U = F^0 \\) の簡略形（エネルギー・ポテンシャル） | \\( 10^{-12} \\) |
| 冪等な \\( n \\) で簡略形の \\( E\_U = 0 \\)；全体のスピン回転で簡略形・sFLL・sAMF が不変 | \\( 10^{-12} \\) |
| dual の和則 \\( \sum\mathrm{Tr}\\,n = N\_{\mathrm e} \\)、\\( \sum H^U\rho = \sum v n \\)（dual・full、共線・非共線） | \\( 10^{-12} \\) |
| \\( \rho \\) 固定の \\( \partial E\_U/\partial\mathbf R \\) の公式と中心差分（dual・full、共線・非共線） | \\( 10^{-8} \\) |
| 拘束：\\( V\lVert T-T\_0\rVert\_F^2 \\) の閉じた形、\\( v = \mathbf h\cdot\boldsymbol\sigma \\) と中心差分（on・on2） | \\( 10^{-8} \\) |
| スピン Zeeman：\\( v = \tfrac12\mathbf b\cdot\boldsymbol\sigma \\) と中心差分 | \\( 10^{-8} \\) |
| 軌道分極の占有の規則が電子数を保存 | \\( 10^{-12} \\) |
| Yukawa の展開と \\( \lambda\to 0 \\) の極限 | \\( 10^{-12} \\)、\\( 10^{-3} \\) |

### 追加の単体テスト〔設計〕

- 2 原子の小さな系で \\( \rho \\)、\\( S \\) を与え、onsite・full・dual の \\( n \\) を手計算と比較（\\( 10^{-14} \\)）。
- 軌道の Zeeman 項：p シェルだけの原子で \\( \Lambda^z \\) が \\( \langle x\vert L\_z\vert y\rangle = -i \\) を再現し、\\( L\_z \\) 固有状態 \\( (x + iy)/\sqrt2 \\) の占有で \\( L\_z = \pm 1 \\) の符号が決まること。
- 拘束 `on` で \\( \mathbf m \parallel \hat{\mathbf e}\_0 \\) のとき \\( E\_{\mathrm{cs}} = 0 \\)、\\( \mathbf h = 0 \\)；反平行のとき \\( E = 2V\lVert\mathbf m\rVert^2 \\)。
- Yukawa 推定：既知の動径関数（例：水素様 3d）で \\( \lambda \\) の二分法が \\( F^0 = U \\) を 0.0005 eV 以内で満たすこと。

### 既存実装との照合（fixture）

| 系 | 設定 | 比べる量 | 許容誤差 |
|---|---|---|---|
| NiO（AFM、共線、マニュアル 40 節の例） | Type 1、dual、Ni 1d に U | `Uhub`、`Utot`、Ni のスピンモーメント、`DM_onsite` の固有値 | 全エネルギー \\( 10^{-6} \\) Hartree、モーメント \\( 10^{-4} \\) μB、占有 \\( 10^{-5} \\) |
| 同上 | onsite・full | `Uhub`、力（有限差分とも比較） | 同上、力 \\( 10^{-5} \\) Hartree/bohr |
| 同上 | Type 2、sFLL・sAMF・cFLL・cAMF、`scf.Yukawa` on/off | `Uhub`、出力の \\( \lambda \\)、\\( J \\)、\\( F^k \\) | \\( 10^{-6} \\) Hartree、\\( \lambda \\) は \\( 10^{-4} \\) |
| FeO（非共線＋スピン軌道、`work/FeO_NC.dat` 相当） | Type 1、dual | `Uhub`、スピン・軌道モーメントとその角度 | \\( 10^{-6} \\) Hartree、\\( 10^{-4} \\) μB、角度 \\( 10^{-3} \\) 度 |
| 同上 | 拘束 on・on2（V = 0.5 eV）、Zeeman スピン・軌道（100 T） | `Ucs`、`Uzs`、`Uzo`、角度 | \\( 10^{-6} \\) Hartree |
| FePt（マニュアル 42 節） | c2n、θ = 0°・90° | `Uele` の差（MAE） | \\( 10^{-6} \\) Hartree |

- 照合は、同じ密度行列から 1 回だけ Hamiltonian を作る単発の比較（SCF の経路の差を除く）と、SCF 収束後の比較の 2 段で行います。
- 軌道分極の促進を使う fixture は、照合モード（4 節の挙動を再現）で比べます。

## 〔未確定〕

- 拘束のポテンシャルで既存実装が使う目標 \\( T\_0 \\) の微分が、正確な勾配と一致するか（不一致なら SCF 収束後の `Ucs` と角度に差が出る）。
- `cFLL`・`cAMF` の XC の平均化が GGA でどう働くか（LDA は確認済み）、`Uxc` のエネルギー評価にも平均化した密度を使うか。
- 共線の dual・onsite での応力の寄与。
- `<Hubbard.U.values` ブロックがないときの U の値（0 になるはずだが未確認）。
- 第二変分法の再開（`c2n`）で、格子密度以外（密度行列、Hamiltonian、占有行列）をどう扱うか。
- 軌道分極の促進の、非共線での軌道添字の転置が意図したものか（照合モードで再現するが、既定の挙動にするかは要相談）。
- 本家 OpenMX 3.9／4.0 と参照した GPU 版で、ここに書いた規則（特に 3.5 節・4 節）が同じか。
