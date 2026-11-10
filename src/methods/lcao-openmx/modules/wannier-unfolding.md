# Wannier 関数・バンドアンフォールディング・局在軌道解析

> 状態：初版（レビュー待ち）。MLWF・CWF・アンフォールディング・軌道最適化の式は公開資料をこのノートの記法で書き直したもの、
> 既定値・帯の選び方・正規化などは既存実装の挙動として確認したものです。
> 式の性質は `tools/check_wannier_unfolding.py` で数値的に確認しています。

## この章の範囲

SCF で得た LCAO の固有状態（[k 点と Bloch 和](kpoints-bloch.md)、[一般化固有値問題](generalized-eigen.md)）と
実空間の行列 \\( H(\mathbf L) \\)、\\( S(\mathbf L) \\)（[疎ブロック形式](../design/sparse-block-format.md)）を入力とする
後処理・解析のうち、次の 5 つを扱います。

| 節 | 機能 | 主な入力キーワード | 公開資料 |
|---|---|---|---|
| 1 | 最大局在 Wannier 関数（MLWF） | `Wannier.*`、`Wannier90.fileout` | マニュアル 47.2、Marzari–Vanderbilt PRB **56**, 12847、Souza–Marzari–Vanderbilt PRB **65**, 035109、Weng–Ozaki–Terakura PRB **79**, 235118 |
| 2 | 最近接 Wannier 関数（CWF） | `CWF.*` | マニュアル 47.1、技術ノート CWF（Ozaki PRB **110**, 125115） |
| 3 | バンドアンフォールディング | `Unfolding.*` | マニュアル 55、Lee–Yamada-Takamura–Ozaki JPCM **25**, 345501（arXiv:1211.4690）、技術ノート unfolding_otw（Lee ほか PRB **102**, 075143） |
| 4 | 自然原子軌道（NAO）・自然電子密度解析と O(N) 軌道解析 | `NBO.*`、`NAO.*` | マニュアル 36、Ohwaki–Otani–Ozaki JCP **140**, 244105 |
| 5 | 基底の軌道最適化（縮約 PAO の生成） | `orbitalOpt.*`、`CntOrb.fileout` | マニュアル 26、技術ノート 5（Ozaki PRB **67**, 155108） |

ここでの「マニュアル」は OpenMX 4.0 のマニュアル（https://www.openmx-square.org/openmx_man4.0/index.html）の節番号です。
SCF そのもの、バンド分散の出力形式（`Band.dispersion`）、cube ファイルの形式は扱いません。

## 記法〔設計〕

[記法・単位・用語](../../../conventions/notation.md) と [k 点と Bloch 和](kpoints-bloch.md) の記法を使います。追加分は次のとおりです。

| 記号 | 意味 |
|---|---|
| \\( N\_k \\) | k 点の数（Born–von Kármán 周期のセル数と同じ） |
| \\( c\_{I\alpha,\mu}(\mathbf k) \\) | LCAO 係数。\\( C(\mathbf k)^\dagger S(\mathbf k) C(\mathbf k) = I \\) |
| \\( \varepsilon\_{\mu\mathbf k} \\) | Kohn–Sham 固有値。\\( \mu\_{\mathrm{F}} \\) は化学ポテンシャル（出力上の「Fermi level」） |
| \\( \psi\_{\mu\mathbf k} \\) | \\( \psi\_{\mu\mathbf k}(\mathbf r) = \sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L}\sum\_{I\alpha} c\_{I\alpha,\mu}(\mathbf k)\\, \phi\_{I\alpha}(\mathbf r - \mathbf R\_I - \mathbf L) \\) |
| \\( w\_{n\mathbf L} \\) | Wannier 関数（セル \\( \mathbf L \\) の \\( n \\) 番目）。\\( N\_W \\) はセルあたりの数 |
| \\( U(\mathbf k) \\) | 帯の空間から Wannier 関数への変換（\\( N\_{\mathrm{band}}\times N\_W \\)、\\( U^\dagger U = I \\)） |

\\( \psi\_{\mu\mathbf k} \\) は 1 セルあたりで規格化されています：\\( \langle\psi\_{\mu\mathbf k}\vert\psi\_{\nu\mathbf k'}\rangle = N\_k\delta\_{\mathbf k\mathbf k'}\delta\_{\mu\nu} \\)。
Bloch 和の位相は **格子ベクトルだけ** に付けます（原子位置を含めない。[k 点と Bloch 和](kpoints-bloch.md)）。
この規約では \\( c(\mathbf k + \mathbf G) = c(\mathbf k) \\) と選べるので、k 点を逆格子ベクトルだけずらしても係数を付け替える必要はありません。

Wannier 関数と、それによる Hamiltonian の行列要素（ホッピング）は、MLWF・CWF に共通して

$$
w\_{n\mathbf L}(\mathbf r) = \frac{1}{N\_k}\sum\_{\mathbf k} e^{-i\mathbf k\cdot\mathbf L}\sum\_{\mu}\psi\_{\mu\mathbf k}(\mathbf r)\\, U\_{\mu n}(\mathbf k),
\qquad
t\_{mn}(\mathbf L) \equiv \langle w\_{m\mathbf 0}\vert\hat H\vert w\_{n\mathbf L}\rangle
{}= \frac{1}{N\_k}\sum\_{\mathbf k} e^{-i\mathbf k\cdot\mathbf L}\left[U(\mathbf k)^\dagger E(\mathbf k) U(\mathbf k)\right]\_{mn}
$$

です（\\( E(\mathbf k) = \mathrm{diag}(\varepsilon\_{\mu\mathbf k}) \\)）。Wannier 補間は \\( H^{W}(\mathbf k) = \sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L}\\, t(\mathbf L) \\) で、
SCF の Bloch 和と同じ符号です。

---

## 1. 最大局在 Wannier 関数（MLWF）

### 1.1 拡がり汎関数〔公開〕

Marzari–Vanderbilt（PRB **56**, 12847）の拡がり

$$
\Omega = \sum\_{n=1}^{N\_W}\left[\langle r^2\rangle\_n - \lvert\langle\mathbf r\rangle\_n\rvert^2\right]
$$

を、k 点格子上の有限差分で評価します。各 k 点で隣接する k 点への差ベクトル \\( \mathbf b \\) と重み \\( w\_b \\) を

$$
\sum\_{\mathbf b} w\_b\\, b\_\alpha b\_\beta = \delta\_{\alpha\beta}
$$

を満たすように選び、重なり行列

$$
M\_{mn}(\mathbf k,\mathbf b) = \langle u\_{m\mathbf k}\vert u\_{n,\mathbf k+\mathbf b}\rangle
{}= \frac{1}{N\_k}\langle\psi\_{m\mathbf k}\vert e^{-i\mathbf b\cdot\mathbf r}\vert\psi\_{n,\mathbf k+\mathbf b}\rangle
$$

（\\( u = e^{-i\mathbf k\cdot\mathbf r}\psi \\)、\\( U(\mathbf k) \\) を掛けた後の表現）から

$$
\bar{\mathbf r}\_n = -\frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\\, \mathbf b\\, \mathrm{Im}\ln M\_{nn}(\mathbf k,\mathbf b),
$$

$$
\Omega\_{\mathrm I} = \frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\left(N\_W - \sum\_{mn}\lvert M\_{mn}\rvert^2\right),\quad
\Omega\_{\mathrm{OD}} = \frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\sum\_{m\ne n}\lvert M\_{mn}\rvert^2,\quad
\Omega\_{\mathrm D} = \frac{1}{N\_k}\sum\_{\mathbf k,\mathbf b} w\_b\sum\_n\left(-\mathrm{Im}\ln M\_{nn} - \mathbf b\cdot\bar{\mathbf r}\_n\right)^2
$$

とし、\\( \Omega = \Omega\_{\mathrm I} + \Omega\_{\mathrm D} + \Omega\_{\mathrm{OD}} \\) です（同論文 Eq. 31–36 に相当）。
\\( \Omega\_{\mathrm I} \\) は \\( U(\mathbf k) \\) の選び方（ゲージ）に依らず、帯の部分空間だけで決まります。

ゲージ依存部分の最小化は、反 Hermite 行列 \\( \Delta W(\mathbf k) \\) で \\( U(\mathbf k) \leftarrow U(\mathbf k)\exp\Delta W(\mathbf k) \\) と更新します。勾配は

$$
G(\mathbf k) = 4\sum\_{\mathbf b} w\_b\left(\mathcal A[R] - \mathcal S[T]\right),\qquad
R\_{mn} = M\_{mn}M\_{nn}^{\*},\quad T\_{mn} = \frac{M\_{mn}}{M\_{nn}}\left(\mathrm{Im}\ln M\_{nn} + \mathbf b\cdot\bar{\mathbf r}\_n\right)
$$

（\\( \mathcal A[X] = (X - X^\dagger)/2 \\)、\\( \mathcal S[X] = (X + X^\dagger)/2i \\)）で、最急降下では
\\( \Delta W = \dfrac{\epsilon}{4 w\_{\mathrm{tot}}}\\, G \\)、\\( w\_{\mathrm{tot}} = \sum\_{\mathbf b} w\_b \\) とします（同論文 Sec. IV）。

### 1.2 帯の絡み合いの解消（disentanglement）〔公開〕

Souza–Marzari–Vanderbilt（PRB **65**, 035109）に従い、外側のエネルギー窓に入る帯から各 k 点で \\( N\_W \\) 次元の部分空間を選び、
\\( \Omega\_{\mathrm I} \\) を最小化します。内側の窓に入る状態は部分空間に必ず含めます（凍結）。
各反復で、隣接 k 点の部分空間への射影から作る行列 \\( Z(\mathbf k) \\) の大きい固有値に対応する \\( N\_W \\) 個のベクトルを選び、
前の反復の \\( Z \\) との線形混合（混合係数 `Wannier.Dis.Mixing.Para`）で振動を抑えます。
収束判定は、隣り合う反復の \\( \Omega\_{\mathrm I} \\) の差（Å²）に `Wannier.Dis.Conv.Criterion` を適用します（マニュアル 47.2.3）。
その後、1.1 節のゲージ依存部分の最小化を行います。

### 1.3 LCAO での行列要素〔公開〕〔規約〕

〔公開〕Weng–Ozaki–Terakura（PRB **79**, 235118）は、MLWF の構成に必要な \\( M \\) と初期射影 \\( A \\) を LCAO 係数と原子軌道間の積分で書きます。
このノートの規約では

$$
M\_{mn}(\mathbf k,\mathbf b) = \sum\_{I\alpha,J\beta}\sum\_{\mathbf L} c\_{I\alpha,m}^{\*}(\mathbf k)\\, c\_{J\beta,n}(\mathbf k + \mathbf b)\\,
e^{i(\mathbf k+\mathbf b)\cdot\mathbf L}\\, e^{-i\mathbf b\cdot\mathbf R\_I}\\,
\left\langle\phi\_{I\alpha}(\mathbf r - \mathbf R\_I)\middle\vert e^{-i\mathbf b\cdot(\mathbf r - \mathbf R\_I)}\middle\vert\phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L)\right\rangle
$$

です（\\( \mathbf k + \mathbf b \\) が k 点格子の外に出ても、上の規約では \\( c(\mathbf k+\mathbf b) = c(\mathbf k + \mathbf b - \mathbf G) \\) をそのまま使えます）。

〔規約〕既存実装は次のように評価します。

- 原子対の積分 \\( \langle\phi\_{I\alpha}\vert e^{-i\mathbf b\cdot(\mathbf r - \mathbf R\_I)}\vert\phi\_{J\beta}\rangle \\) は、SCF の実空間グリッド上の数値積分です
  （[実空間グリッド](grid-poisson-xc.md) の格子点和に体積要素を掛けたもの）。
- 位相の基準を bra 側の原子に置いた評価と、同じ積分を ket 側の対として読み替えた評価の平均を \\( M \\) とします（連続積分なら両者は等しく、グリッド誤差だけが異なる）。
- 非共線では、上の式のスピノル 2 成分の和をとります。

〔公開〕初期射影は、Bloch 状態と、位置 \\( \mathbf p\_n \\) に置いた局在関数（射影子）\\( g\_n \\) の重なりです：

$$
A\_{\mu n}(\mathbf k) = \frac{1}{N\_k}\langle\psi\_{\mu\mathbf k}\vert g\_n\rangle
{}= \sum\_{I\alpha}\sum\_{\mathbf L} e^{-i\mathbf k\cdot\mathbf L}\\, c\_{I\alpha,\mu}^{\*}(\mathbf k)\\, \langle\phi\_{I\alpha}(\mathbf r - \mathbf R\_I - \mathbf L)\vert g\_n(\mathbf r - \mathbf p\_n)\rangle .
$$

射影子 \\( g\_n \\) は PAO（またはその混成）で、重なりは [二中心積分](two-center-overlap-kinetic.md) と同じ方法で計算できます。
帯が孤立している場合の初期 \\( U \\) は \\( A \\) の極分解（\\( A = W\Sigma V^\dagger \\) のとき \\( U = WV^\dagger \\)）で、絡み合う場合は同じ射影で初期部分空間を選びます。

### 1.4 射影子の指定〔公開〕〔規約〕

〔公開〕射影子用の原子種を `Definition.of.Atomic.Species` に通常の原子種と同じ書式で追加し（擬ポテンシャルの欄は書式上必要なだけで使われない）、
`<Wannier.Initial.Projectors` ブロックの各行に「`原子種名-軌道名`、中心の 3 座標、局所 z 軸の 3 成分、局所 x 軸の 3 成分」を書きます（マニュアル 47.2.1）。
各角運動量の動径関数は 1 本だけ使うので、`s1p1d1f1` のような指定を推奨しています。

軌道名と射影子の数は次のとおりです。混成は局所座標系の PAO（[実球面調和関数](pao-basis.md) の並び）の 1 次結合です。

| 軌道名 | 数 | 内容 |
|---|---|---|
| `s`、`px`、`py`、`pz`、`dz2`、`dx2-y2`、`dxy`、`dxz`、`dyz`、`fz3`、`fxz2`、`fyz2`、`fzx2`、`fxyz`、`fx3-3xy2`、`f3yx2-y3` | 1 | その PAO 1 個 |
| `p`、`d`、`f` | 3、5、7 | その \\( l \\) の全成分 |
| `sp`、`sp2`、`sp3` | 2、3、4 | s と p の標準的な混成（\\( sp^3 \\) は各係数の絶対値が \\( 1/2 \\)） |
| `sp3dz2`、`sp3deg` | 5、6 | s、p と d（\\( d\_{z^2} \\)、\\( d\_{x^2-y^2} \\)）を含む混成 |

注：マニュアル 47.2.1 の表では \\( sp^3 \\) 混成の係数が \\( 1/\sqrt 2 \\) と書かれていますが、4 成分の和なので規格化には \\( 1/2 \\) が必要です。
〔規約〕既存実装も \\( 1/2 \\) を使っています。また同じ表の `dxz`、`dyz` の説明が「\\( d\_{xy} \\)」となっているのは誤記で、それぞれ \\( d\_{xz} \\)、\\( d\_{yz} \\) です。

〔規約〕

- 局所座標系は、与えた z 軸・x 軸からオイラー角を求め、\\( l = 1, 2, 3 \\) の実球面調和関数の回転行列で射影子を回します（s は回転不変）。
- 座標の単位 `Wannier.Initial.Projectors.Unit` は `Ang`（既定）、`AU`、`FRAC`。`FRAC` のときは各成分を \\( [0, 1] \\) に折り返してから Cartesian 座標に直します。
- 射影子の総数が `Wannier.Func.Num` と一致しない場合、初期推定を使う（`Wannier.Initial.Guess on`）ならエラーです。

### 1.5 帯の選択〔規約〕

- 窓の端は化学ポテンシャル基準の eV で与えます。各 k 点・各スピンで、外側の窓 \\( [\mu\_{\mathrm F} + E\_{\mathrm{ob}},\ \mu\_{\mathrm F} + E\_{\mathrm{ot}}) \\)（下端を含み上端を含まない）に入る帯を選びます。内側の窓も同じ扱いです。
- 帯の数は k 点ごとに違ってよく、ファイル出力上の帯の数 \\( N\_{\mathrm{win}} \\) はその最大値です。
- ある k 点で外側の窓の帯の数が \\( N\_W \\) 未満、または内側の窓の帯の数が \\( N\_W \\) を超える場合はエラーです。
- 内側の窓の上端と下端の差が \\( 10^{-6} \\) eV 未満なら、内側の窓は使いません（凍結なし）。上端 \\( < \\) 下端はエラー、内側の窓が外側の窓からはみ出すのもエラーです。

### 1.6 b ベクトルの殻〔公開〕〔規約〕

〔公開〕k 点格子の近い順に殻（同じ長さの \\( \mathbf b \\) の組）を探し、\\( \sum w\_b b\_\alpha b\_\beta = \delta\_{\alpha\beta} \\) の 6 本の式を満たす重みを求めます
（Mostofi ほか、Wannier90 の論文 arXiv:0708.0650 Eq. 26 の方法）。探す殻の上限が `Wannier.MaxShells` です（マニュアル 47.2.1）。

〔規約〕殻を 1 つずつ加え、6 × 殻数の係数行列の特異値に \\( 10^{-5} \\) 未満のもの（既存の殻と線形従属）が出たらその殻を捨てて次に進み、
最小二乗解が 6 本の式を各成分 \\( 10^{-5} \\) 以内で満たした時点で終わります。

### 1.7 最小化の手順〔公開〕〔規約〕

〔公開〕`Wannier.Minimizing.Scheme` は 0 = 最急降下（SD）、1 = 共役勾配（CG、割線法による直線探索）、2 = SD の後に CG。
SD の刻みは `Wannier.Minimizing.StepLength`（1.1 節の \\( \epsilon \\)）、CG の割線法の初期刻みと最大回数は
`Wannier.Minimizing.Secant.StepLength`、`Wannier.Minimizing.Secant.Steps` です。
収束判定は勾配のノルムに `Wannier.Minimizing.Conv.Criterion` を適用します。方式 2 では SD と CG がそれぞれ `Wannier.Minimizing.Max.Steps` 回まで回ります（マニュアル 47.2.3）。

〔規約〕\\( \mathrm{Im}\ln M\_{nn} \\) は分枝を前の反復から連続になるように選びます（\\( \lvert M\_{nn}\rvert < 10^{-8} \\) なら分枝の補正だけ）。

### 1.8 キーワードと既定値〔公開〕〔規約〕

既定値は既存実装で確認した値です。マニュアルと異なるものは注記しました。

| キーワード | 既定値 | 意味 |
|---|---|---|
| `Wannier.Func.Calc` | off | MLWF の計算 |
| `Wannier.Func.Num` | 0（実質必須） | \\( N\_W \\) |
| `Wannier.Outer.Window.Bottom` / `.Top` | −10.0 / 10.0 eV | 外側の窓。マニュアルは「既定値なし」 |
| `Wannier.Inner.Window.Bottom` / `.Top` | 0.0 / 0.0 eV | 内側の窓（同じ値なら使わない） |
| `Wannier.Initial.Guess` | on | 射影による初期推定 |
| `Wannier.Initial.Projectors.Unit` | `Ang` | `Ang`、`AU`、`FRAC` |
| `Wannier.Kgrid` | 4 4 4 | MLWF 用の k 点格子。マニュアルは「既定値なし」 |
| `Wannier.MaxShells` | 12 | 殻の探索の上限 |
| `Wannier.Dis.SCF.Max.Steps` | 200 | 絡み合い解消の最大反復 |
| `Wannier.Dis.Conv.Criterion` | \\( 10^{-8} \\) | 同、収束判定（Å²） |
| `Wannier.Dis.Mixing.Para` | 0.5 | 同、混合係数 |
| `Wannier.Minimizing.Scheme` | 0 | 0 = SD、1 = CG、2 = SD→CG |
| `Wannier.Minimizing.StepLength` | 2.0 | SD の刻み |
| `Wannier.Minimizing.Secant.Steps` | 5 | 割線法の最大回数 |
| `Wannier.Minimizing.Secant.StepLength` | 2.0 | 割線法の初期刻み |
| `Wannier.Minimizing.Conv.Criterion` | \\( 10^{-8} \\) | 勾配ノルムの収束判定 |
| `Wannier.Minimizing.Max.Steps` | 200 | 最大反復 |
| `Wannier.Readin.Overlap.Matrix` | off | `.mmn` を読み直して再開。マニュアルは「既定 on」と書くが実装は off |
| `Wannier.Interpolated.Bands` | off | 補間バンド（`Band.dispersion on` が必要。ないとエラー） |
| `Wannier.Function.Plot` | off | MLWF の cube 出力 |
| `Wannier.Function.Plot.SuperCells` | 0 0 0 | 描画するセルの範囲。\\( 2l\_i + 1 \\) が `Wannier.Kgrid` の各成分を超えるとエラー |
| `Wannier90.fileout` | off | Wannier90 用ファイルの出力 |

〔規約〕再開時は外側の窓と k 点格子が保存時と同じであることを確認します。内側の窓・初期推定・収束条件は変えてよいです。

### 1.9 出力〔公開〕〔規約〕

〔公開〕（マニュアル 47.2.2・47.2.5・47.2.7）

| ファイル | 内容 |
|---|---|
| `System.Name.mmn` | \\( M\_{mn}(\mathbf k,\mathbf b) \\)。Wannier90 の形式。2 行目に \\( N\_{\mathrm{win}} \\)、k 点数、\\( \mathbf b \\) の数、スピン数。ブロックはスピン→k 点→\\( \mathbf b \\) の順、各ブロック先頭に k 点番号、\\( \mathbf k+\mathbf b \\) の番号と格子の 3 整数 |
| `System.Name.amn` | \\( A\_{\mu n}(\mathbf k) \\)。Wannier90 の形式（帯番号、Wannier 番号、k 点番号、実部、虚部） |
| `System.Name.eigen` | 化学ポテンシャル、帯数、固有値と固有ベクトル（再開用） |
| `System.Name.eig` | Wannier90 用の固有値（化学ポテンシャル基準の eV） |
| `System.Name.HWR` | \\( t\_{mn}(\mathbf L) \\)。Wigner–Seitz 超セル内の \\( \mathbf L \\) と縮退度、格子ベクトル（bohr）、化学ポテンシャルと値は Hartree |
| `System.Name.Wannier_Band` | 補間バンド（`.Band` と同じ形式） |
| `System.Name.mlwf{s}_{n}_{r,i}.cube` | MLWF の実部・虚部。\\( s \\) はスピン、\\( n \\) は MLWF の番号 |
| `System.Name.win` | Wannier90 の制御ファイル（`Wannier90.fileout on` のとき） |

標準出力には、絡み合い解消（行末 `DISE`、\\( \Omega\_{\mathrm I} \\) とその変化）、最小化（`CONV`、勾配ノルムと \\( \Delta\Omega \\)）、
拡がりの内訳（`SPRD`、\\( \Omega\_{\mathrm I} \\)、\\( \Omega\_{\mathrm D} \\)、\\( \Omega\_{\mathrm{OD}} \\)、合計）、中心（`CENT`、Å と Å²）が出ます。

注：マニュアル 47.2.2 は cube ファイル名を `.mlwfl_…` と書いていますが、〔規約〕実際の名前は `.mlwf{s}_{n}_r.cube` です。
同じく補間バンドの拡張子は `.Wannier_Band` です。

---

## 2. 最近接 Wannier 関数（CWF）

### 2.1 方法〔公開〕

技術ノート CWF（Ozaki PRB **110**, 125115、Eq. 2–18）を、このノートの記法で書き直します。
案内関数（guiding function）\\( Q\_{\mathbf 0 p} \\)（\\( p = 1, \dots, N\_{\mathrm{CWF}} \\)、ホームセルに局在）と窓関数 \\( w(\varepsilon) \\) から射影行列

$$
a\_{\mu p}(\mathbf k) = w(\varepsilon\_{\mu\mathbf k})\\, \frac{1}{N\_k}\langle\psi\_{\mu\mathbf k}\vert Q\_{\mathbf 0 p}\rangle
$$

（\\( N\_{\mathrm{band}}\times N\_{\mathrm{CWF}} \\)、\\( N\_{\mathrm{CWF}} \le N\_{\mathrm{band}} \\)）を作ります。
窓で重み付けた射影関数 \\( L\_{\mathbf R p} \\) と、正規直交な Wannier 関数 \\( W\_{\mathbf R p} \\)（係数 \\( B(\mathbf k) \\)、\\( B^\dagger B = I \\)）の距離の 2 乗和

$$
F[B] = \frac{1}{N\_k}\sum\_{\mathbf k}\left\lVert A(\mathbf k) - B(\mathbf k)\right\rVert\_{\mathrm F}^2
$$

は、\\( A(\mathbf k) \\) の極分解 \\( A = UP \\)（\\( A = W\Sigma V^\dagger \\) の特異値分解から \\( U = WV^\dagger \\)、\\( P = V\Sigma V^\dagger \\)）で \\( B = U \\) としたとき最小で、

$$
F\_{\min} = \frac{1}{N\_k}\sum\_{\mathbf k}\sum\_p\left(\sigma\_{p\mathbf k} - 1\right)^2
$$

です。反復は不要です。性質は次のとおりです。

- \\( A \\) の特異値がすべて正なら \\( U \\) は一意で、Bloch 状態の位相や特異ベクトルの位相・縮退の自由度に依りません（ゲージの問題がない）。
- \\( A^\dagger A \\) が正則なら \\( U = A(A^\dagger A)^{-1/2} \\)、つまり窓付きの Löwdin 直交化と同じです。極分解は \\( (A^\dagger A)^{-1/2} \\) を作らないので、悪条件でも安定です。
- 窓が狭く、ある k 点で窓に入る状態が \\( N\_{\mathrm{CWF}} \\) 個未満になると特異値が 0 になり、\\( U \\) が一意でなくなって対称性が壊れます。窓関数に小さな定数 \\( \delta \\) を足してこれを防ぎます。
- ホッピングと Wannier 補間は「記法」節の式に \\( U \\) を入れたものです。
- \\( N\_{\mathrm{band}} \\) 個の状態から \\( N\_{\mathrm{CWF}} \\) 個を作るので、帯の絡み合いの解消は常に窓関数で自動的に行われます。

### 2.2 窓関数〔公開〕〔規約〕

〔公開〕マニュアル 47.1.1 の 2 種類です。\\( x = \varepsilon - \mu\_{\mathrm F} \\)（化学ポテンシャル基準）とします。

**多項式窓**（`Poly`）：`CWF.disentangling.Erange` に \\( \epsilon\_0 < \epsilon\_1 \le \epsilon\_2 < \epsilon\_3 \\)（eV）を与え、

$$
w\_P(x) = \begin{cases}
0 & x \le \epsilon\_0 \\\\
P\\!\left(\dfrac{\epsilon\_1 - x}{\epsilon\_1 - \epsilon\_0}\right) & \epsilon\_0 < x \le \epsilon\_1 \\\\
1 & \epsilon\_1 < x \le \epsilon\_2 \\\\
P\\!\left(\dfrac{x - \epsilon\_2}{\epsilon\_3 - \epsilon\_2}\right) & \epsilon\_2 < x \le \epsilon\_3 \\\\
0 & \epsilon\_3 < x
\end{cases},\qquad P(y) = 1 - 6y^2 + 8y^3 - 3y^4 .
$$

\\( P(0) = 1 \\)、\\( P'(0) = 0 \\)、\\( P(1) = P'(1) = P''(1) = 0 \\) なので、窓は内側の端で \\( C^1 \\)、外側の端で \\( C^2 \\) です。
外側の区間で厳密に 0 なので、その状態を計算から除けます。〔規約〕多項式窓には \\( \delta \\) を足しません。

**Fermi 窓**（`Fermi`）：`CWF.disentangling.Erange` に \\( \epsilon\_0, \epsilon\_1 \\)、`CWF.disentangling.smearing.kBT` に \\( k\_BT\_0, k\_BT\_1 \\)（eV）を与え、

$$
w\_F(x) = \frac{1}{\left(1 + e^{(\epsilon\_0 - x)/k\_BT\_0}\right)\left(1 + e^{(x - \epsilon\_1)/k\_BT\_1}\right)} + \delta,
\qquad \delta = \texttt{CWF.disentangling.smearing.bound} .
$$

注：技術ノート CWF の Eq. 4 は 2 つの Fermi 関数の和から 1 を引いた形 \\( f\_0 + f\_1 - 1 \\) で、上の積 \\( f\_0 f\_1 \\) とは異なります
（窓の幅が \\( k\_BT \\) より十分広ければほぼ同じ）。〔規約〕既存実装はマニュアルどおり積の形です。

### 2.3 案内関数〔公開〕〔規約〕

`CWF.Guiding.Orbitals` で 3 種類から選びます（マニュアル 47.1.1–47.1.4、技術ノート Sec. III B）。

- **AO**（原子軌道）：PAO そのもの。`<CWF.Guiding.AOs` に原子種ごとに 1 行、基底の \\( (l, \zeta) \\) の順（[PAO 基底](pao-basis.md) の並び）に
  「`1s`、`2s`、…、`1p`、…」の名前と 0/1 を書きます。基底の全 \\( (l,\zeta) \\) を 0 か 1 で明示する必要があります。
  原子あたりの CWF の数は \\( \sum\_{(l,\zeta)\ \mathrm{selected}}(2l+1) \\) です。
- **HO**（混成軌道）：選んだ AO の空間で、その原子のオンサイトの密度行列ブロック \\( \rho\_{\mathbf 0 I\alpha,\mathbf 0 I\beta} \\) を対角化し、
  占有の大きい順に、AO の数と同じ個数の固有ベクトルを使います。選んだ AO が同じなら、AO と HO の CWF は同じ部分空間を張ります。
  〔規約〕密度行列はスピンごとの SCF の密度行列（Fermi 分布の占有）で、窓関数は使いません。PAO はオンサイトで正規直交なので（[PAO 基底](pao-basis.md)）、
  固有ベクトルの Euclid 内積での正規直交性がそのまま関数の正規直交性になります。スピン偏極ではスピンごとに作ります。
- **MO**（埋め込み分子軌道）：`<CWF.MO.Grouped.Atoms` で全原子を群 \\( g \\) に分け、群のブロック
  \\( \Lambda\_g = \sum\_{\mathbf R g'}\rho\_{\mathbf 0 g,\mathbf R g'}\\, S\_{\mathbf R g',\mathbf 0 g} \\)（\\( \rho S \\) の群の対角ブロック）の右特異ベクトルを使います。
  技術ノートの Eq. 31–33 のとおり、\\( \Lambda\_g^\dagger\Lambda\_g \\) を対角化して実の特異ベクトル \\( y\_{g\nu} \\) を得て、
  群内のオンサイト重なり \\( S\_{\mathbf 0 g,\mathbf 0 g} \\) で規格化し、期待値 \\( \langle\bar y\vert\hat H\vert\bar y\rangle \\) を計算します。
  `<CWF.MO.Selection` で群ごとに使う個数を指定します。

〔規約〕MO 案内関数の細部：

- \\( \rho \\) は窓関数 \\( w(\varepsilon) \\) を占有の代わりに使った密度行列です（マニュアル 47.1.4 の記述どおり）。スピン偏極では両スピンの平均（\\( \rho \\) と \\( H \\) とも）を使います。
- 選ぶのは **特異値の大きい順** に先頭から指定個数です。技術ノートは期待値で並べ替えて選ぶと書いていますが、既存実装は特異値の順です。
- `System.Name.CWF_MOs_info` に、群ごとに番号、期待値（化学ポテンシャル基準の eV）、特異値（「population」と表示）を特異値の大きい順に書きます。
- 群番号は 1 から連続していなければなりません（群の数は使われた番号の数で、番号をそのまま添字に使う）。
  `CWF.MO.Selection` の行も群 1 から順に並べます。

〔設計〕案内関数 \\( Q \\) は PAO の 1 次結合なので、\\( \langle\psi\_{\mu\mathbf k}\vert Q\rangle/N\_k = \left[C(\mathbf k)^\dagger S(\mathbf k)\\, q\right]\_{\mu} \\)（\\( q \\) は \\( Q \\) の PAO 係数を全体の並びに置いたもの）で計算できます。

### 2.4 帯の範囲と k 点〔規約〕

- 多項式窓では、窓が 0 でない区間に入る帯の番号の範囲を **SCF の k 点格子の最初の k 点（スピン 1）の固有値** で決め、全 k 点で同じ番号の範囲を使います。
  範囲の帯の数が \\( N\_{\mathrm{CWF}} \\) 未満ならエラーです。
- Fermi 窓では全部の帯（\\( N\_{\mathrm{band}} \\) = 基底関数の数）を使います。
- 周期系（`scf.EigenvalueSolver Band`）では SCF の k 点格子（`scf.Kgrid`）上で計算し、時間反転で対になる k 点は片方だけ計算して重み 2 にします。
  ホッピングは \\( \mathbf L = \sum\_i l\_i\mathbf a\_i \\)、\\( l\_i = -(N\_i - 1)/2, \dots, (N\_i - 1)/2 \\) の範囲（Wigner–Seitz 超セルや縮退度は使わない）で作るため、
  **各方向の k 点数は奇数** でなければなりません（マニュアル 47.1.7）。共線計算ではホッピングの実部だけを保持します。
- 孤立系（`Cluster`）では Γ 点だけで、\\( \mathbf L = \mathbf 0 \\) だけです。
- `CWF.Plot.SuperCells` の \\( 2l\_i + 1 \\) が k 点数を超えるとエラーです。

### 2.5 距離関数と有効電荷の出力〔公開〕〔規約〕

〔規約〕出力される距離関数（`DM_func`）は \\( \sum\_{\mathbf k} \omega\_{\mathbf k}\sum\_p(\sigma\_{p\mathbf k} - 1)^2 / \sum\_{\mathbf k}\omega\_{\mathbf k} \\)（\\( \omega\_{\mathbf k} \\) は k 点の重み）を、
非スピン偏極では 2 倍したものです。併記される「per (CWF*spin)」はそれを \\( 2N\_{\mathrm{CWF}} \\) で割った値です。
マニュアル 47.1.1 は、この値が 0.7 程度を超えると案内関数の対称性が Bloch 状態と合っていない目安だとしています〔公開〕。

〔公開〕CWF の占有（技術ノート Eq. 34）

$$
n\_p = \frac{1}{N\_k}\sum\_{\mathbf k\mu} f(\varepsilon\_{\mu\mathbf k})\\, \lvert u\_{\mu p}(\mathbf k)\rvert^2
$$

から、原子の有効電荷 \\( \kappa\_I = Z\_I^{\mathrm{val}} - \sum\_{p\in I} n\_p \\)（非スピン偏極では \\( n\_p \\) を 2 倍）が得られます。
`.out` ファイルには、CWF ごと・原子ごとの占有（スピン別、和、差）と全体の和が出力されます（マニュアル 47.1.8）。
\\( f \\) は SCF と同じ電子温度の Fermi 分布です。

### 2.6 キーワードと既定値〔公開〕〔規約〕

| キーワード | 既定値 | 意味 |
|---|---|---|
| `CWF.Calc` | off | CWF の計算（`scf.EigenvalueSolver` が `Cluster` か `Band` のとき） |
| `CWF.Guiding.Orbitals` | `HO` | `AO`、`HO`、`MO` |
| `CWF.window.function` | `Poly` | `Poly`、`Fermi` |
| `CWF.disentangling.Erange` | Poly：−1000 −100 0 3／Fermi：−1000 0（eV） | 窓の端。同じキーワードで個数が異なる |
| `CWF.disentangling.smearing.kBT` | 0.3 0.3（eV） | Fermi 窓の \\( k\_BT\_0, k\_BT\_1 \\) |
| `CWF.disentangling.smearing.bound` | \\( 10^{-10} \\) | Fermi 窓の \\( \delta \\)。マニュアルは \\( 10^{-12} \\) と書く（下の注） |
| `<CWF.Guiding.AOs` | — | AO・HO の選択 |
| `<CWF.MO.Grouped.Atoms` / `<CWF.MO.Selection` | — | MO の群と個数 |
| `CWF.Dis.vs.H` | off | 距離とホッピングのテキスト出力 |
| `CWF.fileout.flag` | off | CWF の cube 出力 |
| `CWF.Plot.SuperCells` | 0 0 0 | 描画するセルの範囲 |
| `CWF.fileout.Num` / `<CWF.fileout.cube` | 0 / — | 描画する原子（AO・HO）または群（MO）の番号 |
| `CWF.Wannier90.fileout` | off | Wannier90 用の `.mmn`、`.amn`、`.eig`、`.win` |

注：〔規約〕`Wannier.Func.Calc on` を同時に指定すると、MLWF のキーワードの読み込みで `CWF.disentangling.smearing.bound` の既定値が \\( 10^{-12} \\) で読み直されます。
CWF だけを使う通常の場合の既定値は \\( 10^{-10} \\) です。

### 2.7 出力〔公開〕

（マニュアル 47.1.1・47.1.4・47.1.6・47.1.9・47.1.10）

| ファイル | 内容 |
|---|---|
| `System.Name.CWF.Band` | Wannier 補間したバンド（`Band.dispersion` の経路） |
| `System.Name.CWF.Dis_vs_H` | スピン、原子 GA、原子 GB、\\( l\_1, l\_2, l\_3 \\)、CWF 番号 \\( i, j \\)、中心間距離（Å）、ホッピング（eV）。非共線ではスピンの列がない |
| `System.Name.CWF.Hop` | ホッピングのバイナリ（補助プログラム用） |
| `System.Name.CWF_MOs_info` | MO 案内関数の一覧（2.3 節） |
| `System.Name.cwf_{s}_{atom}_{p}_r.cube` | 共線の CWF（常に実数） |
| `System.Name.cwf_{atom}_{p}.nccube` | 非共線の CWF（\\( \alpha \\) 実部、\\( \alpha \\) 虚部、\\( \beta \\) 実部、\\( \beta \\) 虚部の順） |
| `System.Name_{s}.amn` など | Wannier90 用（スピン添字付き。非共線は添字なし） |

---

## 3. バンドアンフォールディング

### 3.1 定式化〔公開〕

Lee–Yamada-Takamura–Ozaki（JPCM **25**, 345501、arXiv:1211.4690 の Eq. 1–26）を、このノートの記法で書き直します。

超セル（格子ベクトル \\( \mathbf L \\)、原子 \\( I \\)）の固有状態を、より小さい参照セル（格子ベクトル \\( \mathbf r \\)、体積比 \\( n = V\_{\mathrm{super}}/V\_{\mathrm{ref}} \\)）の
Brillouin ゾーンで表します。準備として、

- 超セルの各原子 \\( I \\) に、参照セル内の「原子の種類」\\( g(I) \\)（`Unfolding.Map` の群番号）を割り当てます。
- 位置 \\( \mathbf x \\) を含む参照セルの格子ベクトルを \\( \mathbf r(\mathbf x) = \sum\_i\lfloor\mathbf b\_i^{\mathrm{ref}}\cdot(\mathbf x - \mathbf o)/2\pi\rfloor\\, \mathbf a\_i^{\mathrm{ref}} \\) とします（\\( \mathbf o \\) は参照セルの原点）。

参照セルの k 点 \\( \mathbf k \\) で超セルの固有値問題を解くと（超セルの Bloch 位相 \\( e^{i\mathbf k\cdot\mathbf L} \\) は、\\( \mathbf k \\) を超セルの BZ に折り返した点と同じ）、
状態 \\( \mu \\) の \\( \mathbf k \\) へのアンフォールド重みは

$$
W\_\mu(\mathbf k) = \sum\_{I\alpha}\sum\_{J\beta,\mathbf L}\ \sum\_{M:\\, g(M) = g(J)}
c^{\*}\_{I\alpha,\mu}(\mathbf k)\\, c\_{M\beta,\mu}(\mathbf k)\\,
e^{i\mathbf k\cdot\left[\mathbf r(\mathbf R\_J + \mathbf L) - \mathbf r(\mathbf R\_M)\right]}\\, S\_{I\alpha,J\beta}(\mathbf L)
$$

で、スペクトル重みは \\( W\_\mu(\mathbf k)/n \\) です。和は、ホームセルの軌道 \\( I\alpha \\) と重なりを持つすべての軌道 \\( (J\beta, \mathbf L) \\) について、
それを「同じ種類 \\( g(J) \\) の原子 \\( M \\) の、参照セル \\( \mathbf r(\mathbf R\_J + \mathbf L) \\) への平行移動」と読み替えてとります
（\\( \beta \\) は \\( J \\) と \\( M \\) で同じ軌道番号。同じ群の原子は同じ原子種・同じ基底である必要があります）。
原文の 2 つの近似（参照セルの軌道の位置を超セルの実際の位置で代用する）により、\\( S \\) は超セル計算の重なりそのものです。

性質（`tools/check_wannier_unfolding.py` で確認）：

- 和則：\\( \mathbf k \\) を超セルの逆格子ベクトルだけずらした \\( n \\) 個の点についての和 \\( \sum\_{\mathbf G} W\_\mu(\mathbf K + \mathbf G)/n = 1 \\)（\\( \mathbf c^\dagger S(\mathbf K)\mathbf c = 1 \\) から）。欠陥があっても厳密に成り立ちます。
- 完全な超セル（欠陥なし）では \\( W\_\mu/n \\) は 0 か 1 で、重み 1 の状態のエネルギーは参照セルのバンドに一致します。このとき原点 \\( \mathbf o \\) の選び方に依りません。
- 軌道分解：\\( W\_\mu = \sum\_{M\beta} W\_{\mu,M\beta} \\)（上の式で \\( M\beta \\) を固定した部分和）で、どの軌道が帯を作るかを見られます（原文 Eq. 26）。

注：原文の Eq. 25（および技術ノート unfolding_otw の Eq. 16）は、複素共役を読み替えた側の係数 \\( C\_M \\) に付けています。
このノートの Bloch 規約（\\( C^\dagger S C = 1 \\)、位相は ket 側）で字義どおりに使うと、欠陥のある系で上の和則が成り立ちません（検算で約 1.5 のずれ）。
共役は bra 側（ホームセルの軌道 \\( I\alpha \\)）の係数に付けます。〔規約〕既存実装もそうなっています。

### 3.2 参照セルの原点と対応付け〔公開〕〔規約〕

〔公開〕`<Unfolding.ReferenceVectors` で参照セルを与え（単位は `Atoms.UnitVectors.Unit`、省略すると超セル自身）、
`<Unfolding.Map` で原子ごとに群番号を与えます（正の整数なら番号の値や順序は任意）。原点は `<Unfolding.ReferenceOrigin` で明示できます（マニュアル 55.2・55.3）。

〔規約〕

- `Unfolding.Map` を省略すると、全原子が別々の群になります。群番号は出現順に詰め直されます。
- 原点の既定値：超セルの中心から参照セルの中心を引いた位置に置き、さらに Cartesian の各成分を −0.1 bohr ずらします。
  ある原子（の周期像）の参照セルに対する分数座標が、3 成分とも整数から \\( \tfrac12(0.1,0.1,0.1)\cdot\mathbf b\_i^{\mathrm{ref}}/2\pi \\) 以内にあるときは、もう一度 −0.1 bohr ずらします。
  原点を明示したときはずらしません。
- 同じ原子の異なる周期像が同じ参照セルに入る場合（参照セルが超セルより大きい方向がある場合など）はエラーです。
- 注：マニュアル 55.3 は、原点を「同じ群の原子がなるべく 1 つずつ各参照セルに入る」「原子を含む参照セルの数が最小になる」ように自動決定すると説明していますが、
  確認した版の実装は上の固定規則です。表面・強く歪んだ系では、原点を明示することを勧めます。

### 3.3 k 点・エネルギー・出力〔公開〕〔規約〕

〔公開〕（マニュアル 55.1）

- `<Unfolding.kpoint` に `Unfolding.Nkpoint` 個の点を、**参照セル** の逆格子ベクトルに対する分数座標で並べます（名前、3 座標）。隣り合う点を結んだ折れ線が経路です。
- `Unfolding.desired_totalnkpt` は経路全体の点数の目安です。
- エネルギーの範囲は `Unfolding.LowerBound`、`Unfolding.UpperBound`（化学ポテンシャル基準の eV）。

〔規約〕

- 経路の総長を \\( D \\)（bohr⁻¹）、目安の点数を \\( N\_d \\) として、\\( N\_d \\) が区間数以下なら各区間 2 点（両端）、
  そうでなければ区間 \\( i \\) の点数を \\( \max(1, \lfloor D\_i/\Delta\rfloor) + 1 \\)、\\( \Delta = D/(N\_d - 1) \\) とし、両端を含めて等分します。
  区間のつなぎ目の点は 2 回計算され、横軸の距離は同じ値になります。
- 別形式として `UnfoldBand.Nkpath` と `<UnfoldBand.kpath`（各行：始点 3 座標、終点 3 座標、名前 2 つ）も受け付けます（マニュアルに記載なし）。
- \\( \mu\_{\mathrm F} + E\_{\mathrm{lower}} \le \varepsilon \le \mu\_{\mathrm F} + E\_{\mathrm{upper}} \\)（両端を含む）の状態だけを出力します。
- 固有値問題は \\( S \\) の負の固有値を \\( 10^{-14} \\) に置き換える規則（[一般化固有値問題](generalized-eigen.md) の規則 A′）で解きます。
- 全重みとして \\( \lvert\mathrm{Re}\\, W\_\mu\rvert/n \\) を、軌道分解として \\( \mathrm{Re}\\, W\_{\mu,M\beta}/n \\)（負なら 0 に置き換え）を書きます。
- `scf.EigenvalueSolver` が `Cluster` か `Band` 以外のときは無効になります。軌道最適化を使った計算では縮約後の基底の \\( H \\)、\\( S \\) を使います。

出力ファイル〔公開〕：

| ファイル | 内容 |
|---|---|
| `System.Name.unfold_totup`（`_totdn`、非共線は `_tot`） | k の累積距離（bohr⁻¹）、エネルギー（eV、化学ポテンシャル基準）、全重み |
| `System.Name.unfold_orbup`（`_orbdn`、非共線は `_orb`） | 同じ 2 列の後に軌道ごとの重み（並びは `.out` に記載） |
| `System.Name.unfold_plotexample` | gnuplot の例 |

〔公開〕強度図は重み \\( w \\) を \\( L(k,E) = w/\left[(k/\Delta\_k)^2 + (E/\Delta\_E)^2 + 1\right] \\) でぼかして作ります（マニュアル 55.4、補助プログラム）。
絶対値に物理的な意味はありません。非共線・スピン軌道・DFT+U でも追加のキーワードなしに使えます（マニュアル 55.5）。

### 3.4 光学遷移重みのアンフォールディング〔公開〕

技術ノート unfolding_otw（Lee ほか PRB **102**, 075143）は、Kubo–Greenwood の光学伝導度を k 点と状態ごとの寄与（部分光学伝導度）に分け、
それを 3.1 節の重みで参照セルの BZ に移す方法を述べています。このノートの記法で、超セルの部分光学伝導度

$$
\sigma\_{ab}(\mathbf K, J, \omega) = \frac{-i}{\Omega}\sum\_{J'}\frac{f\_{\mathbf KJ} - f\_{\mathbf KJ'}}{\varepsilon\_{\mathbf KJ} - \varepsilon\_{\mathbf KJ'}}\\,
\frac{\langle\mathbf KJ\vert\hat P\_a\vert\mathbf KJ'\rangle\langle\mathbf KJ'\vert\hat P\_b\vert\mathbf KJ\rangle}{\varepsilon\_{\mathbf KJ} - \varepsilon\_{\mathbf KJ'} + \omega + i\eta}
$$

（\\( \Omega \\) は超セルの体積、縮退時は分数を占有のエネルギー微分で置き換える）に対し、参照セルでの部分光学伝導度は、
スペクトル関数テンソルを 3.1 節と同じ重み \\( W\_J(\mathbf k)/n \\) で移し、体積を参照セルのもの \\( \Omega\_{\mathrm{ref}} = \Omega/n \\) に替えて

$$
\sigma\_{ab}(\mathbf k, j, \omega) = \frac{\Omega}{n\\,\Omega\_{\mathrm{ref}}}\sum\_{\mathbf G}\delta\_{\mathbf k - \mathbf G,\mathbf K}\\, W\_J(\mathbf k)\\, \sigma\_{ab}(\mathbf K, J, \omega)
{}= \sum\_{\mathbf G}\delta\_{\mathbf k - \mathbf G,\mathbf K}\\, W\_J(\mathbf k)\\, \sigma\_{ab}(\mathbf K, J, \omega)
$$

となります（原文 Eq. 15–17 をこのノートの \\( W \\) の正規化で書き直したもの。原文の係数 \\( (L/l)^2 \\) は原文の \\( L \\)、\\( l \\) の定義に依存します）。
周波数区間で積分したものをバンド図に重ねて表示します。
運動量行列要素は LCAO 係数と原子軌道間の \\( \langle\phi\vert\nabla\vert\phi\rangle \\) で書け、軌道 4 つの組に分解できます（原文 Eq. 19–24）。

〔規約〕確認した版には、アンフォールドした光学遷移重みを出力する入力キーワードはありません（光学伝導度とアンフォールディングはそれぞれ独立に計算できる）。

---

## 4. 自然原子軌道（NAO）解析と O(N) 軌道解析

### 4.1 方法〔公開〕

Weinhold らの自然結合軌道法（Reed–Weinstock–Weinhold, JCP **83**, 735 (1985)）の第 1 段階である NAO と自然電子密度（natural population, NP）を計算します（マニュアル 36）。
密度行列 \\( P \\)（スピンごと）と重なり \\( S \\) から、次の順に作ります。

1. 原子 \\( A \\)・角運動量 \\( l \\) ごとに、\\( (2l+1) \\) 個の成分について平均した \\( P\_{Al} \\)、\\( S\_{Al} \\)（\\( \zeta \\) について \\( N\_l\times N\_l \\)）で
   \\( P\_{Al}\mathbf n = w\\, S\_{Al}\mathbf n \\) を解き、前 NAO と占有 \\( w \\) を得ます。
2. 占有がしきい値以上のものを自然最小基底（NMB）、残りを自然 Rydberg 基底（NRB）に分けます。
3. NMB を占有で重み付けた対称直交化（OWSO）で原子間に直交化し、NRB を NMB に対して Schmidt 直交化してから NRB どうしを直交化します。
4. 直交化した基底で原子ブロックの密度行列を再び \\( (A, l) \\) ごとに対角化して NAO を得ます。原子の NP はその対角和です。

〔公開〕Ohwaki–Otani–Ozaki（JCP **140**, 244105）は、O(N) 法の打ち切りクラスターの考え方で、
指定した中心原子のまわりのクラスターの \\( P \\)、\\( S \\) だけから NAO と NP を求め、大規模系の局所的な解析を可能にしました。

### 4.2 入力と規約〔公開〕〔規約〕

| キーワード | 既定値 | 意味 |
|---|---|---|
| `NBO.switch` | `off` | `on1`：Cluster 計算で全系の密度行列を使う。`on2`：Krylov（O(N)）計算で中心原子だけ |
| `NAO.threshold` | 0.85 | NMB と NRB を分ける占有のしきい値（スピン軌道あたりの電子数） |
| `NBO.Num.CenterAtoms` | 1 | `on2` の中心原子の数（全原子数を超えるとエラー） |
| `<NBO.CenterAtoms` | — | 中心原子の番号（`on2` のときだけ有効） |
| `NAO.only` | on | NAO・NP だけ計算 |
| `NHO.fileout`、`NBO.fileout` | off | 自然混成軌道・自然結合軌道の cube 出力 |
| `NBO.SmallCell.Switch`、`<NBO.SmallCell` | off | 出力範囲を分数座標で制限 |

〔規約〕

- `on1` は `scf.EigenvalueSolver Cluster`、`on2` は `Krylov` のときだけ動きます。組み合わせが違うとメッセージを出して何もしません。
  `on3`、`on4`（周期系用）は受け付けますが、確認した版では処理が無効化されています。
- しきい値は、スピン偏極（共線）のとき `NAO.threshold`/2 として、スピンごとの占有と比べます（非スピン偏極では `NAO.threshold` のまま）。
- 原子ブロックの \\( S\_{Al} \\) の固有値に \\( 10^{-14} \\) 未満が出ると停止します。
- `on2` のクラスターは、Krylov 法と同じ打ち切りクラスター（中心原子、基底の重なる原子、さらにその外側の原子）です（[未確定](#未確定) を参照）。
- 結果は標準出力に、原子ごとの NP とその総和、原子ごとの NAO の占有・エネルギー・PAO 展開係数として書かれます。

---

## 5. 基底の軌道最適化（縮約 PAO の生成）

### 5.1 方法〔公開〕

技術ノート 5（Ozaki PRB **67**, 155108、Eq. 1–12）とマニュアル 26 の方法です。基底関数を、原始 PAO（同じ原子・同じ \\( l \\) の動径関数）の 1 次結合で表し、
その係数（縮約係数）\\( a \\) を変分的に最適化します。このノートの記法で

$$
\phi\_{I,l\zeta m}(\mathbf r) = \sum\_{p=0}^{N^{\mathrm{prim}}\_l - 1} a\_{I,l\zeta,p}\\, \chi\_{I,lpm}(\mathbf r),\qquad
\chi\_{I,lpm}(\mathbf r) = R^{(I)}\_{lp}(r)\\, S\_{lm}(\hat{\mathbf r})
$$

とします。\\( a \\) は \\( m \\) に依らない（制限付き最適化）ので、全エネルギーは系の回転に対して不変です。
LCAO 係数 \\( c \\) で SCF を解いた後、力の定理と同じ考え方で

$$
\frac{\partial E\_{\mathrm{tot}}}{\partial a\_{I\alpha,p}} = 2\sum\_{\sigma}\sum\_{J\beta}
\left(\rho^{\sigma}\_{I\alpha,J\beta}\\, \langle\chi\_{Ip}\vert\hat H^{\sigma}\vert\phi\_{J\beta}\rangle - E^{\sigma}\_{I\alpha,J\beta}\\, \langle\chi\_{Ip}\vert\phi\_{J\beta}\rangle\right)
$$

が得られます（\\( \rho \\) は密度行列、\\( E \\) はエネルギー密度行列、いずれも縮約後の基底で、実空間の全原子対・周期像について和をとる。非スピン偏極ではスピンの和を 2 倍とする）。
この勾配は縮約後の基底の \\( \rho \\)、\\( E \\) と、原始基底の \\( H \\)、\\( S \\) だけで決まり、追加の微分は要りません。
制限付き最適化では \\( (2l+1) \\) 個の \\( m \\) についての勾配の和を、その \\( (l,\zeta) \\) のすべての \\( m \\) の勾配とします。

〔検算〕固定した Hamiltonian の有限クラスターで、この勾配がバンドエネルギーの差分商と \\( 10^{-8} \\) 程度で一致することを確認しました（`tools/check_wannier_unfolding.py`）。

### 5.2 指定〔公開〕〔規約〕

〔公開〕基底を `H5.0-s4>1`、`C5.0-s4>1p4>1` のように「原始数 `>` 縮約数」で書きます（[PAO 基底](pao-basis.md)）。
原始関数は PAO ファイルの先頭から \\( N^{\mathrm{prim}}\_l \\) 本です。`orbitalOpt.Method` で、原子ごと（`Atoms`）か原子種ごと（`Species`、同じ原子種名の原子で係数を共有）かを選びます。

〔規約〕

- `orbitalOpt.Method` は `Off`（既定）、`Atoms`、`Species`、`Atoms2`、`Species2`。`Atoms2`・`Species2` は、各 \\( l \\) の **最後の** \\( \zeta \\) の縮約係数だけを最適化し、残りを初期値に固定します（マニュアルに記載なし）。
- `Species` では、原子種の勾配をその原子種の全原子の勾配の和とします。
- 縮約係数は更新のたびに \\( \sum\_{pq}a\_p a\_q\\, S^{(I)}\_{pq} = 1 \\)（オンサイトの原始関数の重なり）となるように規格化します。PAO は正規直交なので実質 \\( \sum\_p a\_p^2 = 1 \\) です。
- 軌道最適化の 1 ステップは、`orbitalOpt.scf.maxIter` 回の SCF 反復ごと、またはそれより前に SCF が収束したとき（エネルギー変化が `scf.criterion` 未満かつ残差ノルムが \\( 10^{-7} \\) 未満）に行います。
- 収束判定の「勾配ノルム」は、勾配の各成分の 2 乗和（平方根をとらない。\\( m \\) 成分ごとに数える）です。これが `orbitalOpt.criterion` 未満になるか、ステップ数が `orbitalOpt.Opt.maxIter` を超えたら、最後にもう 1 回 `orbitalOpt.scf.maxIter` 回の SCF を回して終わります。
- MD・構造最適化では、`orbitalOpt.per.MDIter` ステップごと（1, 1+P, 1+2P, … ステップ目）に最適化をやり直します。
- 併用できない機能：NEGF（入力エラー）、`RMM-DIISH`・`RMM-DIISV`（[密度ミキシング](density-mixing.md)）、応力（[全エネルギーと力](energy-force.md)）。

### 5.3 更新の規則〔公開〕〔規約〕

〔公開〕`orbitalOpt.StartPulay` より前のステップは最急降下、以降は準 Newton 法（`EF`：固有ベクトル追跡、`DIIS`）で、`orbitalOpt.HistoryPulay` 個の履歴を使います（マニュアル 26）。

〔規約〕履歴の数 \\( d = \min(\texttt{HistoryPulay},\ \mathrm{step} - \texttt{StartPulay} + 1) \\)（1 未満なら 1）として、

| モード | \\( d = 1 \\)（最急降下） | \\( d > 1 \\) |
|---|---|---|
| `Atoms` | \\( a \leftarrow a - \lambda g \\)、\\( \lambda = \min(4\\, s/\lVert g\rVert^2,\ 0.05) \\) | DIIS の外挿 \\( \bar a - \gamma\bar g \\)。\\( \gamma \\) は勾配ノルムが \\( 10^{-2}, 10^{-3}, 10^{-4}, 10^{-5} \\) を超えるとき 0.1、0.2、0.4、0.5、それ以下で 0.7（`EF` の指定でも DIIS） |
| `Species` | \\( a \leftarrow a - (s/\max\_i\lvert g\_i\rvert)\\, g \\) | まず DIIS で \\( \bar a \\)、\\( \bar g \\) を作り、`EF` なら BFGS で更新した Hessian（固有値を 0.1 以上に制限）で準 Newton ステップを求め、成分の最大値が 0.001 を超えれば 0.001 に縮めて \\( \bar a \\) から引く。`DIIS` なら \\( \bar a - \gamma\bar g \\)（\\( \gamma \\) は 0.05、0.10、0.15、0.20、0.40） |

ここで \\( s \\) は `orbitalOpt.SD.step`、\\( \lVert g\rVert^2 \\) は上の勾配ノルムです。

### 5.4 初期係数〔公開〕〔規約〕

〔公開〕技術ノート 5 の Eq. 11–12：原子ごとに近接原子を含むクラスターの Kohn–Sham 方程式を非自己無撞着に解き、Fermi 分布で重み付けて初期係数を作ります。

〔規約〕`orbitalOpt.InitCoes` は `Symmetrical`（既定）、`Free`、`Simple`。

- `Symmetrical`・`Free`（確認した版では同じ処理）：原子とその近接原子のクラスターを原始基底で対角化し、クラスターの電子数から決めた局所化学ポテンシャルと
  電子温度 1500 K の Fermi 分布で重み付けて、その原子の \\( l \\) ごとのオンサイト密度行列（\\( m \\) について和）を作ります。
  対角化した固有ベクトルを占有の大きい順に \\( \zeta = 0, 1, \dots \\) の縮約係数とし、絶対値最大の成分が正になるよう符号をそろえます
  （縮退を解くため、対角に \\( -10^{-9}(N\_l - p)^2 \\) を加える）。
- `Simple`：\\( a\_{l\zeta,p} = \delta\_{\zeta p} \\)（\\( \zeta \\) 番目の原始関数そのもの）。

### 5.5 キーワードと既定値〔公開〕〔規約〕

| キーワード | 既定値 | 意味 |
|---|---|---|
| `orbitalOpt.Method` | `Off` | `Off`、`Atoms`、`Species`、`Atoms2`、`Species2` |
| `orbitalOpt.Opt.Method` | `DIIS` | `DIIS`、`EF` |
| `orbitalOpt.InitCoes` | `Symmetrical` | 初期係数 |
| `orbitalOpt.SD.step` | 0.001 | 最急降下の刻み \\( s \\) |
| `orbitalOpt.HistoryPulay` | 20 | 履歴の数。マニュアルは 15 と書く |
| `orbitalOpt.StartPulay` | 10 | 準 Newton を始めるステップ。マニュアルは 1 と書く |
| `orbitalOpt.scf.maxIter` | 40 | 1 ステップあたりの SCF 反復 |
| `orbitalOpt.Opt.maxIter` | 100 | 最大ステップ数 |
| `orbitalOpt.per.MDIter` | 1000000 | MD でやり直す間隔 |
| `orbitalOpt.criterion` | \\( 10^{-4} \\) | 勾配ノルムのしきい値（(Hartree/bohr)² と表示される） |
| `CntOrb.fileout` | off | 縮約 PAO の出力（`Atoms` か `Species` のときだけ） |
| `Num.CntOrb.Atoms` / `<Atoms.Cont.Orbitals` | 1 / — | 出力する原子の数と番号 |

### 5.6 縮約 PAO ファイルの出力〔公開〕〔規約〕

〔公開〕`CntOrb.fileout on` で、指定した原子ごとに `{原子種名}_{原子番号}.pao` を出力し、そのまま基底として使えます（マニュアル 26）。

〔規約〕ファイルの構成は次のとおりです。

1. 説明の行（元の PAO ファイル名、基底の指定、入力ファイル名）。
2. `<Contraction.coefficients` … `Contraction.coefficients>`：各行 `Atom= L= Mul= p= 係数`（\\( m = 0 \\) の成分だけ。原始数を超える \\( p \\) は 0）。
3. 元の PAO ファイルの `<pseudo.atomic.orbitals.L=0` の直前までの写し。
4. 基底に使う \\( l \\) の動径関数の表：先頭の \\( N^{\mathrm{cont}}\_l \\) 列が縮約後の \\( \sum\_p a\_{l\zeta,p}R\_{lp}(r) \\)、
   残りの列（\\( \zeta \ge N^{\mathrm{cont}}\_l \\)）は元の原始関数 \\( R\_{l\zeta} \\) そのもの（縮約後の関数に対して直交化しない）。
5. 基底に使わない \\( l \\)（\\( l\_{\max} \\) を超える \\( l \\)）の表は元のまま。

読み込み側（[PAO 基底](pao-basis.md)）は表の先頭から \\( N\_l \\) 本を使うので、縮約数以下の \\( N\_l \\) で指定すれば最適化した基底になります。

---

## 6. 端のケース〔規約〕〔設計〕

- **窓に状態が足りない**：MLWF で外側の窓の帯数が \\( N\_W \\) 未満の k 点がある、CWF の多項式窓の帯の範囲が \\( N\_{\mathrm{CWF}} \\) 未満、はいずれもエラー。
  CWF で Fermi 窓が狭く \\( \delta \\) が小さいと、特異値が \\( \delta \\) 程度になり、CWF は対称性を保つものの補間の質が落ちます（技術ノート Sec. IV A の Si の例）。
- **案内関数が帯の対称性と合わない**：\\( F \\) が大きくなり、補間バンドが乱れます。2 つの温度 \\( k\_BT\_0 \ne k\_BT\_1 \\) で、捨てたい帯を急峻に、含めたい帯を緩やかに扱えます（技術ノート Fig. 6）。
- **HO で占有のほぼ 0 の AO を選ぶ**：オンサイト密度行列の固有値が 0 付近で縮退し、選んでいない AO の方向と混ざりうる。〔設計〕再実装では、選んだ AO の部分空間に制限した行列を対角化します。
- **CWF の k 点数が偶数**：ホッピングの \\( \mathbf L \\) の範囲が対称にならないため、既存実装は奇数だけを受け付けます。
- **アンフォールディングで群に異なる原子種を混ぜる**：軌道数・軌道の並びが合わず、式が定義されません。〔設計〕入力検査でエラーにします。不純物は別の群にします（マニュアル 55.2）。
- **参照セルの境界上の原子**：原点のわずかな違いで \\( \mathbf r(\cdot) \\) が変わり、欠陥系では重みが変わります。完全結晶では原点に依りません。
- **非直交基底での負の軌道重み**：軌道分解の重みは負になりうる（全重みと和則は正しい）。既存実装は出力時だけ 0 に置き換えます。
- **軌道最適化の縮約数が原始数に等しい**：`s1>1` は `s1` と同じ。最適化する自由度は規格化を除いて \\( N^{\mathrm{prim}} - 1 \\) です。

## 7. 設計：モジュール化した再実装の要件〔設計〕

- **共通の入力**：k 点ごとの \\( \varepsilon\_{\mu\mathbf k} \\)、\\( C(\mathbf k) \\)、\\( S(\mathbf k) \\)、化学ポテンシャル、実空間の \\( S(\mathbf L) \\)、\\( H(\mathbf L) \\)、密度行列。
  SCF の内部状態に依存せず、[疎ブロック形式](../design/sparse-block-format.md) と k 点の結果から再計算できるようにします。
- **Wannier の核**：(a) 射影行列 \\( A(\mathbf k) \\) を作る部分（案内関数・射影子・窓関数を差し替え可能）、(b) 極分解（CWF）、
  (c) \\( M(\mathbf k,\mathbf b) \\) と MV/SMV の最小化（MLWF）、(d) \\( U(\mathbf k) \\) からホッピング・補間・占有を作る部分、に分けます。
  (d) は MLWF と CWF で共通にし、\\( \mathbf L \\) の集合（対称な箱か Wigner–Seitz 超セルと縮退度か）を引数にします。
- **窓関数** は純関数として持ち、多項式窓・Fermi 窓（積の形）を既定とし、技術ノートの \\( f\_0 + f\_1 - 1 \\) 形も選べるようにします。
- **\\( M \\) の積分**：既存実装は実空間グリッドの和ですが、再実装では二中心積分（[二中心積分](two-center-overlap-kinetic.md) の Fourier 法に位相 \\( e^{-i\mathbf b\cdot\mathbf r} \\) を加えたもの）で評価してもよいです。差は fixture の許容誤差で吸収します。
- **ファイル形式**：`.mmn`、`.amn`、`.eig` は Wannier90 の形式に合わせ、Wannier90 で同じ結果が出ることをテストにします。
- **アンフォールディング** は「参照セル・群・原点 → 原子ごとの \\( \mathbf r(\cdot) \\) と群」の前処理と、「k ごとの重み」の計算に分けます。
  重みは疎ブロックの \\( S(\mathbf L) \\) を原子対で走査すれば O(N) で、状態ごとに独立に並列化できます。原点は既存実装の規則を既定とし、明示指定を推奨します。
- **NAO** は「原子ブロックの抽出 → 前 NAO → OWSO → NAO」の純粋な線形代数として、入力の \\( P \\)、\\( S \\) を全系またはクラスターから与えられるようにします（O(N) 版はクラスターを作る部分だけが違う）。
- **軌道最適化** は SCF ループの外側の最適化として、勾配計算（\\( \rho \\)、\\( E \\) と原始基底の \\( H \\)、\\( S \\) から）と更新規則（SD、DIIS、BFGS）を分けます。
  縮約は [PAO 基底](pao-basis.md) の動径関数表の 1 次結合として表し、縮約後の表を二中心積分・グリッド計算にそのまま渡します。
- **既定値の扱い**：マニュアルと実装で異なる既定値（`Wannier.Kgrid`、窓、`Wannier.Readin.Overlap.Matrix`、`CWF.disentangling.smearing.bound`、`orbitalOpt.HistoryPulay`、`orbitalOpt.StartPulay`）は、互換モードで実装の値に合わせます。

## 8. テスト〔設計〕

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 多項式窓：\\( P(0)=1 \\)、\\( P'(0)=0 \\)、\\( P(1)=P'(1)=P''(1)=0 \\)、区間の端で連続 | \\( 10^{-12} \\) |
| 解析 | 極分解：\\( F[U] = N\_k^{-1}\sum(\sigma-1)^2 \\)、任意の部分ユニタリー \\( B \\) で \\( F[B] \ge F[U] \\) | \\( 10^{-12} \\) |
| 解析 | \\( A \\) が満階数なら \\( U = A(A^\dagger A)^{-1/2} \\) | \\( 10^{-8} \\) |
| 不変量 | Bloch 状態にランダムな位相を掛けても CWF のホッピングが不変 | \\( 10^{-12} \\) |
| 不変量 | 共線・時間反転対称な系で CWF のホッピングが実数 | \\( 10^{-10} \\) |
| 解析 | 全軌道を案内関数、\\( w \equiv 1 \\)：格子上の k で補間バンドが元のバンドに一致 | \\( 10^{-10} \\) |
| 不変量 | MLWF：\\( \Omega\_{\mathrm I} \\) がゲージ（ランダムなユニタリー）に依らない | \\( 10^{-12} \\) |
| 解析 | MLWF：完全な組（点状の正規直交軌道）では \\( \Omega\_{\mathrm I}=0 \\)、最小化で \\( \Omega\to 0 \\)、中心が軌道の位置に一致 | \\( \Omega<10^{-5} \\)、中心 \\( 10^{-3} \\) |
| 解析 | b ベクトルの重みが \\( \sum w\_b b\_\alpha b\_\beta = \delta\_{\alpha\beta} \\) を満たす（単純立方・fcc・六方） | \\( 10^{-10} \\) |
| 不変量 | アンフォールディングの和則 \\( \sum\_{\mathbf G}W/n = 1 \\)（欠陥・変位あり） | \\( 10^{-10} \\) |
| 解析 | 完全な超セル：重み 0/1、重み 1 のエネルギーが参照セルのバンドに一致、原点に依らない | \\( 10^{-10} \\) |
| 解析 | 参照セル = 超セル（`Unfolding.ReferenceVectors` 省略）：全状態の全重みが 1 | \\( 10^{-10} \\) |
| 解析 | 軌道最適化の勾配が全エネルギー（または固定 \\( H \\) のバンドエネルギー）の差分商と一致 | 相対 \\( 10^{-6} \\) |
| 解析 | 縮約係数の規格化 \\( a^{\mathsf T}S^{(I)}a = 1 \\)、出力 `.pao` の縮約列が \\( \sum\_p a\_pR\_p \\) | \\( 10^{-12} \\) |
| 照合 | Si（ダイヤモンド）MLWF（マニュアルの例、8 個）：\\( \Omega\_{\mathrm I} \\)、合計の拡がり、中心 | 拡がり 1%、中心 \\( 10^{-3} \\) Å |
| 照合 | GaAs の CWF（マニュアル 47.1.1 の条件）：距離関数、ホッピング、CWF の占有 | ホッピング \\( 10^{-4} \\) eV、占有 \\( 10^{-5} \\) |
| 照合 | 2 次元 SiC の 2×2 超セル（完全・Si 空孔）のアンフォールド重み（マニュアル 55.2） | \\( 10^{-6} \\) |
| 照合 | エチレンカーボネートの NP（マニュアル 36 の例） | \\( 10^{-4} \\) 電子 |
| 照合 | メタンの軌道最適化（マニュアル 26 の例）：勾配ノルムの履歴の最初の数ステップ、最終全エネルギー | 勾配ノルム相対 \\( 10^{-3} \\)、エネルギー \\( 10^{-6} \\) Hartree |
| 相互 | MLWF・CWF の `.mmn`、`.amn`、`.eig` を Wannier90 に読ませた拡がり・補間バンドが本体の結果と一致 | 拡がり 1% |

検算：`uv run --with numpy --with scipy python tools/check_wannier_unfolding.py` は、1 次元の Gaussian 軌道模型（非直交）で
アンフォールディングの和則・完全結晶での 0/1 重み・原点非依存性、原文の共役の位置では和則が破れること、
CWF の極分解の最小性・ゲージ不変性・Löwdin との一致・補間、多項式窓の滑らかさ、MLWF の \\( \Omega\_{\mathrm I} \\) の不変性と最小化、
軌道最適化の勾配を確認します（すべて合格）。

## 〔未確定〕

- MLWF の計算が `scf.EigenvalueSolver Cluster`（孤立系）でも使えるか、使える場合の b ベクトルの扱い。
- MLWF の絡み合い解消（\\( Z \\) 行列の作り方、凍結状態の扱い、初期部分空間の選び方）の細部が Wannier90 と一致するか。fixture での確認が必要です。
- `.HWR` の Wigner–Seitz 超セルの作り方（Wannier90 と同じ範囲・縮退度の判定か）。
- MLWF の CG の割線法と、方式 2 の SD→CG の切り替えの細部（ステップの受理条件）。
- CWF の多項式窓で、帯の番号の範囲を 1 つの k 点の固有値で決める既存実装の規則を、再実装の既定として引き継ぐか（k 点ごとに窓の中の帯を全部使う方が自然）。
- CWF の非共線計算での HO・MO 案内関数（複素になりうる）の作り方と、スピン添字の扱い。
- CWF の有効電荷の出力で、Cluster と Band の正規化（k 点の重み）が同じか。
- NAO：スピン偏極でしきい値を半分にする規則の意図（スピンごとの占有の最大値は 1 なので、0.425 は低すぎる可能性がある）。
- NAO の O(N) 版（`on2`）のクラスターの正確な範囲（Krylov 法の章で定める打ち切りクラスターとの関係）と、NHO・NBO（`NAO.only off`）の処理。この章では NAO・NP までを対象とします。
- 軌道最適化：`orbitalOpt.Force.Skip` の意味、`Atoms2`・`Species2` の初期係数の作り方（最後の \\( \zeta \\) に最大占有の固有ベクトルを入れる）の意図、
  `orbitalOpt.criterion` の単位表示（(Hartree/bohr)² と表示されるが、実際には縮約係数についての勾配の 2 乗和）。
- アンフォールディングの原点の自動決定について、マニュアル 55.3 の説明（探索による決定）と確認した版の実装（固定規則）のどちらを再実装の既定にするか。
