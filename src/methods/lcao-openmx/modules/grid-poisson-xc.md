# 実空間グリッド・Poisson・交換相関

> 状態：試作（レビュー待ち）。式は技術ノート 1（Ver. 1.3、Eq. 15, 21, 47–58, 63–80）、技術ノート 2、
> マニュアル 11・14 節をこのノートの記法で書き直し、公開資料に記載のない規則は既存実装の挙動として〔規約〕で補いました。
> 数値の検算は `tools/check_grid_poisson_xc.py` で行っています。

この章は、単位胞を等間隔に分割した実空間格子の作り方、その上での密度の表現、
原子密度からの差 \\( \delta n \\) に対する Poisson 方程式、交換相関（XC）ポテンシャル、
格子上のポテンシャルからの行列要素を扱います。中性原子ポテンシャル \\( V\_{\mathrm{na}} \\) は
[中性原子ポテンシャル](neutral-atom.md)、エネルギーの分け方と力の全体は [全エネルギーと力](energy-force.md) を参照します。

## 1. 実空間格子

### 1.1 格子点〔公開〕

技術ノート 1（Eq. 48–50）に従い、格子ベクトル \\( \mathbf a\_i \\)（\\( i = 1,2,3 \\)）をそれぞれ \\( N\_i \\) 等分します。

$$
\mathbf h\_i = \frac{\mathbf a\_i}{N\_i},\qquad
\Delta V = \lvert \mathbf h\_1\cdot(\mathbf h\_2\times\mathbf h\_3)\rvert = \frac{\Omega}{N\_1N\_2N\_3},\qquad
\mathbf r\_{\mathbf n} = \mathbf O + \sum\_i n\_i \mathbf h\_i,\quad n\_i = 0,\dots,N\_i-1
$$

\\( \Omega \\) は単位胞の体積、\\( \mathbf O \\) は格子の原点です。
\\( \mathbf h\_i \\) の逆格子ベクトルは \\( \mathbf g\_i = N\_i \mathbf b\_i \\)（\\( \mathbf b\_i \\) は単位胞の逆格子ベクトル、
\\( \mathbf a\_i\cdot\mathbf b\_j = 2\pi\delta\_{ij} \\)）です。
格子上の関数は周期 \\( N\_i \\) で巡回するものとして扱い、周期像の位置 \\( \mathbf r\_{\mathbf n} + \mathbf L \\) の値は \\( \mathbf r\_{\mathbf n} \\) の値と同じです。

### 1.2 分割数 \\( N\_i \\) の決め方〔規約〕

`scf.energycutoff` の値 \\( E \\)（Ry、既定値 220）から、各方向について目標の分割数

$$
N\_i^{\ast} = \frac{2\sqrt{E}}{\lVert\mathbf b\_i\rVert}
$$

を求め、\\( 2^p3^q5^r7^s \\)（\\( p,q,r,s \ge 0 \\)）の形の整数のうち \\( \lvert \ln N - \ln N\_i^{\ast}\rvert \\) が最小のものを \\( N\_i \\) とします。
最も近い数を対数尺度で選ぶので、\\( N\_i < N\_i^{\ast} \\) となる（実際のカットオフが指定値より小さくなる）こともあります。

方向ごとの「実際のカットオフ」は

$$
E\_i = \left(\frac{\lVert\mathbf g\_i\rVert}{2}\right)^2 = \left(\frac{N\_i\lVert\mathbf b\_i\rVert}{2}\right)^2 \quad(\mathrm{Ry})
$$

で、\\( \lVert\mathbf g\_i\rVert/2 \\) はその方向で格子が表せる最大の波数（Nyquist 波数）です。
直方体の単位胞では \\( N\_i^{\ast} = \sqrt{E}\\,\lVert\mathbf a\_i\rVert/\pi \\)、格子間隔は約 \\( \pi/\sqrt{E} \\) bohr です
（150 Ry で 0.257 bohr）。例（`tools/check_grid_poisson_xc.py` の出力）：

| 単位胞（bohr） | \\( E \\)（Ry） | \\( (N\_1,N\_2,N\_3) \\) | \\( E\_i \\)（Ry） |
|---|---|---|---|
| 立方 10 | 150 | (40, 40, 40) | 157.91 |
| 立方 10 | 220 | (48, 48, 48) | 227.40 |
| 直方体 5 × 7.3 × 21 | 200 | (24, 32, 96) | 227.40, 189.65, 206.25 |
| fcc（\\( a = 5.4 \\)） | 220 | (15, 15, 15) | 228.46 |

> **技術ノート 1 Eq. 47 についての注記**：同式は方向ごとのカットオフを \\( \tfrac12\lVert\mathbf g\_i\rVert^2 \\) と書いていますが、
> 既存実装が採用し出力する値は上の \\( E\_i = (\lVert\mathbf g\_i\rVert/2)^2 \\) Ry（Hartree では \\( \tfrac18\lVert\mathbf g\_i\rVert^2 \\)）で、4 倍異なります。
> 再実装では \\( N\_i \\) の決め方を上の規則で定義し、Eq. 47 は「三方向のカットオフをそろえる」という趣旨としてだけ参照します。

〔公開〕マニュアル 14.1 節は、多くの場合 150〜200 Ry が適当で、ポテンシャル面が平坦な場合やスピン配置間の小さなエネルギー差では
300 Ry 以上が必要になることがあるとしています。

### 1.3 分割数の固定（`scf.Ngrid`）

〔公開〕マニュアル 14.2 節：`scf.Ngrid n1 n2 n3` を与えると \\( N\_i \\) をその値に固定し、`scf.energycutoff` は使いません。
格子定数を変えながらエネルギー曲線を描くとき、\\( N\_i \\) の変化による不連続を避けるためのものです。

〔規約〕

- `scf.Ngrid` の各値は \\( 2^p3^q5^r7^s \\) の形でなければならず、そうでなければ入力エラーとします。3 つとも 0（既定）のとき固定しません。
- 固定したとき、報告用の要求カットオフは \\( (E\_1+E\_2+E\_3)/3 \\) に置き換えます。
- 格子の最適化（セルの最適化）では、最初に決めた \\( N\_i \\) を以後固定します。

### 1.4 格子の原点〔規約〕

既定では、格子点の重心と原子位置の算術平均が一致するように原点を置きます。

$$
\mathbf O = \frac{1}{N\_{\mathrm{at}}}\sum\_I \mathbf R\_I - \sum\_i \frac{N\_i-1}{2}\\,\mathbf h\_i
$$

- 原点は計算の最初（MD・構造最適化の最初のステップ）に一度だけ決め、以後原子が動いても変えません。
- `scf.fixed.grid x y z`（bohr）を与えると、原点をその値に固定します（〔公開〕マニュアル 14.3 節：
  一連の計算で原子と格子の相対位置を保ち、格子由来の誤差をそろえるため。値は標準出力に表示される原点をそのまま使う）。
- 〔規約〕ESM（`ESM.switch`）のときは、法線方向の成分について原点の決め方が変わります（原子位置の算術平均を格子間隔の整数倍に切り捨ててから上の式に入れる）。
  規則は [ESM・外部電場・内殻準位](esm-field-corehole.md) 2.2 節にあります。

### 1.5 孤立系（クラスター）の箱〔規約〕

入力に `<Atoms.UnitVectors` ブロックがない場合、孤立系として直方体の箱を自動で作ります
（バンド計算、分率座標、ESM ではエラーとします）。

1. 各原子 \\( I \\) について、その種の PAO ファイルに書かれたカットオフ半径（`radial.cutoff.pao`）\\( r\_{c,I} \\) を使い、
   Cartesian 各成分 \\( \alpha \\) で \\( \min\_I (R\_{I\alpha} - r\_{c,I}) \\) と \\( \max\_I (R\_{I\alpha} + r\_{c,I}) \\) を求めます。
2. 辺の長さをその差の 1.10 倍とし、格子ベクトルを Cartesian 軸に平行に取ります。
3. 格子間隔 \\( h = \pi/\sqrt{E + 0.001} \\)（\\( E \\) は Ry）として、\\( \lfloor L/h\rfloor \\) 以上で \\( jh \ge L \\) を満たす
   最小の \\( 2^p3^q5^r7^s \\) 型の整数 \\( j \\) を探し、辺の長さを \\( jh \\) に広げます。

こうして作った箱は 1.2 節の規則でちょうど \\( N\_i = j \\) を与えます（スクリプトで確認）。
箱はそのまま周期境界条件で扱います。Poisson 方程式の特別な扱いは 3.3 節を参照してください。

## 2. 格子上の密度

### 2.1 電子密度〔公開〕

技術ノート 1（Eq. 6, 15）：基底関数 \\( \phi\_{I\alpha}(\mathbf r) = R\_{l\zeta}(r) S\_{lm}(\hat{\mathbf r}) \\)（[PAO 基底](pao-basis.md)）と
密度行列 \\( \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\) から

$$
n\_\sigma(\mathbf r) = \sum\_{I\alpha}\sum\_{J\beta,\mathbf L} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\,
\phi\_{I\alpha}(\mathbf r - \mathbf R\_I)\\, \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L)
$$

を格子点 \\( \mathbf r\_{\mathbf n} \\) で直接評価します。\\( n = n\_\uparrow + n\_\downarrow \\)。

〔規約〕

- 各原子の周期像について、\\( \lVert \mathbf r\_{\mathbf n} - \mathbf R\_I - \mathbf L\rVert \le r\_{c,I} \\) を満たす格子点の一覧を作ります
  （\\( r\_{c,I} \\) は基底指定の PAO カットオフ半径）。原子対 \\( (I, J+\mathbf L) \\) の寄与はこの 2 つの球の共通部分の格子点だけで計算します。
- スピン非分極では \\( n\_\uparrow = n\_\downarrow = n/2 \\) として 2 成分を保持します。
- 格子上の密度を電子数で規格化し直すことはしません（\\( \Delta V\sum n \\) と電子数の差は離散化誤差としてそのまま残ります）。

### 2.2 原子密度の重ね合わせと \\( \delta n \\)〔公開〕

技術ノート 1（Eq. 6）：

$$
n^{(a)}(\mathbf r) = \sum\_{I,\mathbf L} n^{(a)}\_I(\lVert\mathbf r - \mathbf R\_I - \mathbf L\rVert),\qquad
\delta n(\mathbf r) = n(\mathbf r) - n^{(a)}(\mathbf r)
$$

〔規約〕

- \\( n^{(a)}\_I \\) は PAO ファイルに収められた原子の価電子密度（動径関数）です。読み込み時に、PAO の対数メッシュ上の和で
  \\( 4\pi\int r^2 n^{(a)}\_I\\,dr \\) を求め、価電子数（擬ポテンシャルの価電荷 \\( Z\_I \\)）に一致するよう定数倍します。
- 格子上では、1 節と同じ「\\( r \le r\_{c,I} \\) の格子点」にだけ値を置きます。
- 初期密度として \\( n^{(a)} \\) を使うときのスピンへの振り分けは [密度行列と電荷ミキシング](density-mixing.md) で扱います。

### 2.3 部分内殻補正（PCC）密度〔公開〕〔規約〕

〔公開〕技術ノート 1（Eq. 15）：XC には \\( n\_{\mathrm{pcc}} \\) を各スピンに半分ずつ加えた密度を使います。

$$
E\_{\mathrm{xc}} = \int d^3r\\, \bigl(n(\mathbf r) + n\_{\mathrm{pcc}}(\mathbf r)\bigr)\\,
\varepsilon\_{\mathrm{xc}}\bigl(n\_\uparrow + \tfrac12 n\_{\mathrm{pcc}},\\, n\_\downarrow + \tfrac12 n\_{\mathrm{pcc}}\bigr)
$$

〔規約〕

- \\( n\_{\mathrm{pcc}} = \sum\_{I,\mathbf L} n^{\mathrm{pcc}}\_I \\) で、\\( n^{\mathrm{pcc}}\_I \\) は擬ポテンシャル（VPS）ファイルの部分内殻密度です。
  格子上では原子密度と同じく \\( r \le r\_{c,I} \\) の格子点にだけ置きます。
- `scf.partialCoreCorrection` は常に on でなければならず、off はエラーです（PCC を持たない擬ポテンシャルでは \\( n^{\mathrm{pcc}}\_I = 0 \\)）。
- スピン分極計算では、`<scf.pcc.opencore` ブロックで種ごとに \\( +1 \\)・\\( -1 \\)・\\( 0 \\) を指定でき、
  \\( +1 \\)（\\( -1 \\)）の種の \\( n^{\mathrm{pcc}}\_I \\) は半分ずつでなく全部を上（下）スピンに加えます（開殻内殻の扱い）。

### 2.4 非共線スピン〔公開〕〔規約〕

〔公開〕技術ノート 2（Eq. 3, 22–25, 34–35）：格子点ごとに 2×2 の密度行列

$$
\underline n(\mathbf r) = \begin{pmatrix} n\_{11} & n\_{12}\\\\ n\_{21} & n\_{22}\end{pmatrix}
{}= \tfrac12\left(n\\,\sigma\_0 + \mathbf m\cdot\boldsymbol\sigma\right),\qquad n\_{21} = \overline{n\_{12}}
$$

を持ち、局所的に対角化して局所スピン座標系の密度 \\( n'\_{\uparrow}, n'\_{\downarrow} \\) と磁化の向き \\( (\theta,\varphi) \\) を得ます。
このノートの記法では

$$
\mathbf m = (2\\,\mathrm{Re}\\,n\_{12},\\; -2\\,\mathrm{Im}\\,n\_{12},\\; n\_{11}-n\_{22}),\qquad
n'\_{\uparrow,\downarrow} = \tfrac12\left(n \pm \lVert\mathbf m\rVert\right),
$$

\\( \theta, \varphi \\) は \\( \hat{\mathbf m} \\) の極角・方位角（\\( \hat{\mathbf m} = (\sin\theta\cos\varphi, \sin\theta\sin\varphi, \cos\theta) \\)）です。

〔規約〕\\( \lVert\mathbf m\rVert < 10^{-14} \\) のときは \\( \theta = \pi/2,\ \varphi = 0 \\) とします（このとき \\( n'\_\uparrow = n'\_\downarrow \\) なので結果に影響しません）。
\\( \theta \in [0,\pi] \\)、\\( \varphi \in [-\pi/2, 3\pi/2) \\) の範囲で求めます。

## 3. Hartree ポテンシャル（FFT による Poisson 方程式）

### 3.1 周期系〔公開〕

技術ノート 1（Eq. 52–54）：\\( \delta n \\) の離散 Fourier 係数

$$
\delta\tilde n(\mathbf G) = \frac{1}{N\_1N\_2N\_3}\sum\_{\mathbf n} \delta n(\mathbf r\_{\mathbf n})\\, e^{-i\mathbf G\cdot\mathbf r\_{\mathbf n}},\qquad
\mathbf G = \sum\_i k\_i \mathbf b\_i
$$

から

$$
\delta\tilde V\_H(\mathbf G) = \frac{4\pi}{\lVert\mathbf G\rVert^2}\\,\delta\tilde n(\mathbf G),\qquad
\delta V\_H(\mathbf r\_{\mathbf n}) = \sum\_{\mathbf G} \delta\tilde V\_H(\mathbf G)\\, e^{i\mathbf G\cdot\mathbf r\_{\mathbf n}}
$$

とします（Hartree 原子単位、\\( \delta V\_H \\) は電子 1 個が感じるポテンシャルエネルギーで、電子密度を正に取る規約）。
\\( \delta\tilde n \\) は離散 Fourier 係数で、このノートの連続 Fourier 変換 \\( \tilde f \\)（[記法](../../../conventions/notation.md)）とは
\\( \delta\tilde n(\mathbf G) \approx \tilde{\delta n}(\mathbf G)/\Omega \\) の関係です。
原点 \\( \mathbf O \\) による位相 \\( e^{-i\mathbf G\cdot\mathbf O} \\) は順変換と逆変換で打ち消すので、指数の \\( \mathbf r\_{\mathbf n} \\) は \\( \sum\_i n\_i\mathbf h\_i \\) に置き換えて構いません。

〔規約〕

- 整数 \\( k\_i \\) は \\( 0 \le \kappa\_i < N\_i \\) の添字 \\( \kappa\_i \\) から、\\( \kappa\_i < \lfloor N\_i/2\rfloor \\) なら \\( k\_i = \kappa\_i \\)、
  それ以外は \\( k\_i = \kappa\_i - N\_i \\) とします。逆変換の結果は実部だけを使います。
  - \\( N\_i \\) が偶数のとき、これは通常の対称な選び方（Nyquist 成分を \\( -N\_i/2 \\) 側に置く）と同じで、結果は実数です。
  - \\( N\_i \\) が奇数のとき、\\( \kappa\_i = (N\_i-1)/2 \\) が \\( k\_i = -(N\_i+1)/2 \\) に写り、対称な選び方
    （\\( k\_i \in [-(N\_i-1)/2, (N\_i-1)/2] \\)）とは最高波数の成分だけが異なります。共役な成分の対に異なる \\( \lVert\mathbf G\rVert \\) が割り当てられるため
    逆変換は虚部を持ち、実部を取ることはその対のカーネルを平均することに当たります。
    差は分解能の限界の成分だけに現れ、カットオフを上げると消えます（例：一辺 10 bohr、\\( N = 15 \\)、\\( e^{-3r^2} \\) で
    \\( E\_{\delta\mathrm H} \\) の差 \\( 5\times10^{-5} \\) Ha、\\( N=16 \\) では 0）。
  - 〔設計〕再実装はこの写し方をそのまま採用し（fixture との照合のため）、対称な選び方を選択肢として持ちます。
- \\( \mathbf G = 0 \\) 成分は 0 にします。すなわち \\( \delta V\_H \\) の単位胞平均は 0 です。
  系が電荷を持つ場合（`scf.system.charge`）、これは一様な補償電荷を置いたことに当たります。
- 前向き変換は \\( e^{-i\cdots} \\)、逆変換は \\( e^{+i\cdots} \\) で、正規化 \\( 1/(N\_1N\_2N\_3) \\) は前向き側に付けます。

### 3.2 Hartree エネルギー〔公開〕

技術ノート 1（Eq. 21, 51）：

$$
E\_{\delta\mathrm{H}} = \frac12\int\_\Omega \delta n\\, \delta V\_H\\, d^3r \approx \frac{\Delta V}{2}\sum\_{\mathbf n} \delta n(\mathbf r\_{\mathbf n})\\, \delta V\_H(\mathbf r\_{\mathbf n})
{}= \frac{\Omega}{2}\sum\_{\mathbf G\neq 0}\frac{4\pi}{\lVert\mathbf G\rVert^2}\lvert\delta\tilde n(\mathbf G)\rvert^2
$$

最後の等号は離散 Parseval の等式で、格子上では厳密です。

### 3.3 孤立系：球形 Coulomb カットオフ〔規約〕

`scf.coulomb.cutoff on`（既定 off）のとき、カーネルを

$$
\delta\tilde V\_H(\mathbf G) = \frac{4\pi}{G^2}\left(1 - \cos(G R\_c)\right)\delta\tilde n(\mathbf G)\quad(\mathbf G\ne 0),\qquad
\delta\tilde V\_H(0) = 2\pi R\_c^2\\,\delta\tilde n(0)
$$

に置き換えます。\\( R\_c \\) は単位胞の最短の格子ベクトルの長さの半分です。これは距離 \\( R\_c \\) で切った Coulomb 核
\\( \Theta(R\_c - r)/r \\) の Fourier 変換で、電荷が半径 \\( R\_c/2 \\) 程度の範囲に収まっていれば、
その範囲での孤立系のポテンシャルを周期像の影響なしに与えます（\\( \mathbf G=0 \\) の値は \\( G\to0 \\) の極限）。
既定（off）では、自動で作った孤立系の箱（1.5 節）でも 3.1 節の周期系の式をそのまま使います。
\\( \delta n \\) は中性原子の差なので双極子より高次の寄与が主で、周期像の影響は箱の余白で抑えられる、という扱いです。

〔設計〕ESM 法（`ESM.switch`）による表面系の Poisson 方程式と外部電場は [ESM・外部電場・内殻準位](esm-field-corehole.md) で、NEGF の電極境界条件は [NEGF による電子輸送](negf.md) 7 節で扱い、この章の範囲外とします。

## 4. 交換相関

### 4.1 汎関数の一覧〔公開〕〔規約〕

〔公開〕マニュアル 11 節：`scf.XcType` で選びます。

| `scf.XcType` | 汎関数 | スピン |
|---|---|---|
| `LDA` | Ceperley–Alder（Perdew–Zunger の補間式、PRB **23**, 5048） | `scf.SpinPolarization off` のみ |
| `LSDA-CA` | 同上のスピン分極版 | off・on・NC |
| `LSDA-PW` | Perdew–Wang 92 の相関（PRB **45**, 13244）＋ Slater 交換 | off・on・NC |
| `GGA-PBE` | Perdew–Burke–Ernzerhof（PRL **77**, 3865） | off・on・NC |

〔規約〕`LDA` とスピン分極（on・NC）の組み合わせは入力エラーです。入力は他に `HSE`・`NLX` という値も受け付けますが、
〔規約〕格子上の XC の評価器はこれらに対して交換相関を計算せず（入力は受理されるが実質的に到達不能）、意味のある計算になりません。
再実装は「未対応」として停止します（[厳密交換・分散力補正](exx-vdw.md) の「既存実装の状況」）。LDA+U の一部の二重計数規則では、XC の評価に両スピンの平均密度を使う（スピン分極を 0 とみなす）
変種がありますが、これは DFT+U の章で扱います。

### 4.2 各汎関数の式〔公開〕

上の原論文の式を使います。このノートでの書き方：\\( r\_s = (3/(4\pi n))^{1/3} \\)、\\( \zeta = (n\_\uparrow - n\_\downarrow)/n \\)、

$$
f(\zeta) = \frac{(1+\zeta)^{4/3} + (1-\zeta)^{4/3} - 2}{2^{4/3}-2},\qquad
\varepsilon\_x^{\mathrm{unif}}(n) = -\frac34\left(\frac{3}{\pi}\right)^{1/3} n^{1/3}
$$

- 交換（全汎関数共通）はスピンスケーリング \\( E\_x[n\_\uparrow,n\_\downarrow] = \tfrac12\left(E\_x[2n\_\uparrow] + E\_x[2n\_\downarrow]\right) \\) で、
  LDA では粒子あたり \\( \varepsilon\_x = \sum\_\sigma n\_\sigma\\,\varepsilon\_x^{\mathrm{unif}}(2n\_\sigma)/n \\) です。
  `LDA`・`LSDA-CA` では同じ量を \\( \varepsilon\_x^{\mathrm{P}} + f(\zeta)(\varepsilon\_x^{\mathrm F} - \varepsilon\_x^{\mathrm P}) \\)
  （\\( \varepsilon\_x^{\mathrm F} = 2^{1/3}\varepsilon\_x^{\mathrm P} \\)）と書いても同じです。
- PZ81 相関：非分極（P）・完全分極（F）について、\\( r\_s \ge 1 \\) で \\( \gamma/(1+\beta\_1\sqrt{r\_s}+\beta\_2 r\_s) \\)、
  \\( r\_s < 1 \\) で \\( A\ln r\_s + B + C r\_s\ln r\_s + D r\_s \\)、両者を \\( f(\zeta) \\) で補間します。
  係数は P：\\( (\gamma,\beta\_1,\beta\_2,A,B,C,D) = (-0.1423, 1.0529, 0.3334, 0.0311, -0.048, 0.0020, -0.0116) \\)、
  F：\\( (-0.0843, 1.3981, 0.2611, 0.01555, -0.0269, 0.0007, -0.0048) \\)。
- PW92 相関：原論文の Table I の 3 組の係数（\\( \varepsilon\_c(r\_s,0) \\)、\\( \varepsilon\_c(r\_s,1) \\)、\\( -\alpha\_c(r\_s) \\)）と
  \\( \varepsilon\_c = \varepsilon\_c(r\_s,0) + \alpha\_c f(\zeta)(1-\zeta^4)/f''(0) + [\varepsilon\_c(r\_s,1)-\varepsilon\_c(r\_s,0)] f(\zeta)\zeta^4 \\)。
  〔規約〕\\( f''(0) = 1.70992093416137 \\)（厳密値 \\( 8/(9(2^{4/3}-2)) \\)）を使います（libxc の `LDA_C_PW` は 1.709921）。
- PBE：交換の増強因子 \\( F\_x(s) = 1+\kappa-\kappa/(1+\mu s^2/\kappa) \\)、\\( \kappa = 0.804 \\)、\\( \mu = \beta\pi^2/3 \\)、
  \\( s = \lvert\nabla n\rvert/(2k\_F n) \\)（スピンスケーリングで \\( 2n\_\sigma \\) に適用）。
  相関は PW92 の \\( \varepsilon\_c \\) に勾配補正 \\( H(r\_s,\zeta,t) \\) を加えたもの、\\( \beta = 0.06672455060314922 \\)、
  \\( \gamma = (1-\ln 2)/\pi^2 \\)、\\( t = \lvert\nabla n\rvert/(2\phi k\_s n) \\)、\\( \phi = \tfrac12[(1+\zeta)^{2/3}+(1-\zeta)^{2/3}] \\)。

これらの式は libxc 7.1.2 と相対誤差 \\( 10^{-15} \\)（LDA）、\\( 3\times10^{-7} \\)（PBE。libxc 側の PW92 係数の桁数の差）で一致することを
スクリプトで確認しています。

### 4.3 小さな密度の扱い〔規約〕

物理量に影響しうる閾値として、既存実装は次を使います（密度は bohr⁻³、PCC を加えた後の値）。

| 汎関数 | 規則 |
|---|---|
| `LDA` | \\( n \le 10^{-15} \\) で \\( \varepsilon\_{xc} = v\_{xc} = 0 \\) |
| `LSDA-CA` | \\( n < 10^{-15} \\) で 0。\\( \lvert\zeta\rvert > 1 \\) なら \\( \zeta = \pm(1-10^{-15}) \\) に制限 |
| `LSDA-PW` | \\( n < 10^{-14} \\) で 0。相関では \\( \zeta > 0.99 \\) を 0.99 に制限（\\( \zeta < -1 \\) の場合だけ \\( -0.99 \\) に置換） |
| `GGA-PBE` | 価電子密度 \\( n\_\uparrow + n\_\downarrow < 10^{-14} \\) の点では \\( \varepsilon\_{xc} = v\_{xc} = 0 \\)、勾配も 0。それ以外では各スピン密度を \\( 0.5\times10^{-14} \\) 以上に、\\( \lvert\nabla n\rvert \\) と \\( 2\lvert\nabla n\_\sigma\rvert \\) を \\( 10^{-14} \\) 以上に切り上げ、\\( \zeta \\) は `LSDA-PW` と同じ制限 |

`GGA-PBE` の判定だけは PCC を加える前の価電子密度で行います。〔設計〕再実装では、fixture と比較するときに同じ閾値を使い、
それ以外では、閾値の違いは \\( 10^{-14} \\) bohr⁻³ 程度の密度の点にしか影響しないことをテストで確認すれば十分とします。

### 4.4 格子上の XC エネルギーとポテンシャル〔公開〕

技術ノート 1（Eq. 55, 58）：格子点ごとに \\( \tilde n\_\sigma = n\_\sigma + \tfrac12 n\_{\mathrm{pcc}} \\) を作り

$$
E\_{\mathrm{xc}}^{\mathrm{grid}} = \Delta V\sum\_{\mathbf n}\left(n + n\_{\mathrm{pcc}}\right)\varepsilon\_{\mathrm{xc}}(\tilde n\_\uparrow,\tilde n\_\downarrow,\nabla\tilde n\_\uparrow,\nabla\tilde n\_\downarrow),\qquad
v\_{\mathrm{xc},\sigma} = \frac{\partial e\_{\mathrm{xc}}}{\partial \tilde n\_\sigma} - \nabla\cdot\frac{\partial e\_{\mathrm{xc}}}{\partial \nabla\tilde n\_\sigma}
$$

とします（\\( e\_{\mathrm{xc}} = (n+n\_{\mathrm{pcc}})\varepsilon\_{\mathrm{xc}} \\) は体積あたりのエネルギー）。LDA では第 2 項はありません。

### 4.5 GGA の勾配：27 点の多項式当てはめ〔公開〕〔規約〕

〔公開〕マニュアル 11 節：密度勾配は、中心を含む 27 個の近傍格子点への実空間の多項式当てはめで求めます。

〔規約〕具体的には次のとおりです。

1. 中心の格子点からの変位 \\( \mathbf d\_{jkl} = j\mathbf h\_1 + k\mathbf h\_2 + l\mathbf h\_3 \\)（\\( j,k,l\in\\{-1,0,1\\} \\)）の 27 点を取ります。
2. Cartesian 座標 \\( \mathbf d = (x,y,z) \\) の 27 個の単項式 \\( x^p y^q z^s \\)（\\( p,q,s\in\\{0,1,2\\} \\)）の線形結合を、27 点の値に
   ちょうど一致させます（27×27 の連立方程式。係数行列は格子だけで決まるので一度だけ逆行列を作る）。
3. 勾配は \\( x, y, z \\) の 1 次の係数です。

この方法は全次数 2 以下の多項式に対して厳密で、直交格子では各方向の中心差分 \\( (f\_{+} - f\_{-})/(2h) \\) に一致します（スクリプトで確認）。
勾配を取る密度は PCC を含む \\( \tilde n\_\sigma \\) です。

ポテンシャルの第 2 項 \\( \nabla\cdot(\partial e\_{\mathrm{xc}}/\partial\nabla\tilde n\_\sigma) \\) も、各格子点で求めた
\\( \partial e\_{\mathrm{xc}}/\partial\nabla\tilde n\_\sigma \\) の 3 成分に同じ当てはめを適用して計算します。
したがってポテンシャルは中心から各方向 ±2 点の値に依存します。

〔設計〕FFT による勾配も可能ですが、既存実装との照合のため、既定は上の 27 点法とします。

### 4.6 非共線 XC〔公開〕〔規約〕

〔公開〕技術ノート 2（Eq. 18–19）：局所座標系の密度 \\( n'\_\uparrow, n'\_\downarrow \\) で
\\( v\_{\mathrm{xc},\uparrow}, v\_{\mathrm{xc},\downarrow} \\) を求め、全体座標系に戻します。
\\( \bar v = \tfrac12(v\_{\uparrow}+v\_{\downarrow}) \\)、\\( \Delta v = \tfrac12(v\_\uparrow - v\_\downarrow) \\) として

$$
\underline V\_{\mathrm{xc}} = \bar v\\,\sigma\_0 + \Delta v\\,\hat{\mathbf m}\cdot\boldsymbol\sigma
{}= \begin{pmatrix} \bar v + \Delta v\cos\theta & \Delta v\sin\theta\\, e^{-i\varphi}\\\\ \Delta v\sin\theta\\, e^{i\varphi} & \bar v - \Delta v\cos\theta\end{pmatrix}
$$

です。\\( \mathrm{Tr}(\underline V\_{\mathrm{xc}}\\,\delta\underline n) = v\_\uparrow\delta n'\_\uparrow + v\_\downarrow\delta n'\_\downarrow \\)（1 次）が成り立つことをスクリプトで確認しています。

〔規約〕

- GGA では局所座標系の密度 \\( n'\_\uparrow(\mathbf r), n'\_\downarrow(\mathbf r) \\) を格子上のスカラー場とみなして 4.5 節の勾配を取ります
  （磁化の向きの空間変化は勾配に入れない近似）。
- PCC は \\( n'\_\uparrow, n'\_\downarrow \\) に半分ずつ加えます（開殻内殻の種では 2.3 節の規則で片方へ）。
- Hartree ポテンシャルは対角成分（11 と 22）にだけ加えます。

### 4.7 原子密度の XC エネルギーの精密評価（\\( E\_{\mathrm{xc}}^{(0)} \\) 補正）〔規約〕

格子上の和 \\( E\_{\mathrm{xc}}^{\mathrm{grid}} \\) は原子と格子の相対位置に依存します（エッグボックス効果）。
既存実装は既定で（`scf.Exc0.correction.flag on`）、原子密度の重ね合わせ \\( \hat n = n^{(a)} + n\_{\mathrm{pcc}} \\) の
CA-LDA エネルギーを格子から差し引き、原子中心の細かい球面積分で足し戻します。

$$
E\_{\mathrm{xc}} = \Delta V\sum\_{\mathbf n}\Bigl[(n+n\_{\mathrm{pcc}})\varepsilon\_{\mathrm{xc}}(\tilde n\_\uparrow,\tilde n\_\downarrow,\dots)
{}- \hat n\\,\varepsilon^{\mathrm{CA}}\_{xc}(\hat n)\Bigr] + E\_{\mathrm{xc}}^{(0)}
$$

$$
E\_{\mathrm{xc}}^{(0)} = \sum\_I \int\_{\lVert\mathbf r-\mathbf R\_I\rVert < r\_{c,I}} \hat n\_I(\mathbf r - \mathbf R\_I)\\, \varepsilon^{\mathrm{CA}}\_{xc}\bigl(\hat n(\mathbf r)\bigr)\\, d^3r,\qquad
\hat n\_I = n^{(a)}\_I + n^{\mathrm{pcc}}\_I
$$

- \\( \varepsilon^{\mathrm{CA}}\_{xc} \\) は `scf.XcType` に関係なく非分極の CA-LDA（4.2 節、閾値 \\( 10^{-15} \\)）です。
  \\( \hat n\_I \\) の台がそれぞれの \\( r\_{c,I} \\) に収まるので、\\( \sum\_I \hat n\_I = \hat n \\) により第 2 項は \\( \int\hat n\\,\varepsilon^{\mathrm{CA}}(\hat n) \\) の分割になります。
- 積分は、原子 \\( I \\) を中心に、動径 \\( [0, r\_{c,I}] \\) の 150 点 Gauss–Legendre 則と、角度の Lebedev 則（`scf.Num.Leb.Grid`、既定 590 点）の積で行います。
  \\( \hat n(\mathbf r) \\) は、\\( I \\) と原子間距離 \\( r\_{c,I}+r\_{c,J} \\) 以内の原子 \\( J \\)（周期像を含む）の \\( \hat n\_J \\) の和です。
- \\( \hat n \\) は SCF で変わらないので、ポテンシャル（ハミルトニアン）は変わりません。エネルギーと力だけが変わります。
- off にすると \\( E\_{\mathrm{xc}} = E\_{\mathrm{xc}}^{\mathrm{grid}} \\) です。

〔設計〕再実装は既定を on とし、on・off の両方を fixture と比較します。

## 5. 格子上のポテンシャルの行列要素

### 5.1 有効ポテンシャル〔公開〕〔規約〕

〔公開〕技術ノート 1（Eq. 30–31）：格子で扱う局所ポテンシャルは

$$
V\_\sigma^{\mathrm{grid}}(\mathbf r) = \delta V\_H(\mathbf r) + v\_{\mathrm{xc},\sigma}(\mathbf r)\\;\bigl[+\\, V\_{\mathrm{na}}(\mathbf r)\bigr]
$$

〔規約〕\\( V\_{\mathrm{na}} \\) は `scf.ProExpn.VNA off` のときだけ格子に載せます（既定の on では射影子展開で別に計算）。
外部電場を使う場合はその一様電場のポテンシャルも加えます（[ESM・外部電場・内殻準位](esm-field-corehole.md) 3 節で扱う）。非共線では \\( \delta V\_H \\)（と \\( V\_{\mathrm{na}} \\)）を対角の 11・22 成分に加え、
12 成分は \\( \underline V\_{\mathrm{xc}} \\) だけです。

### 5.2 行列要素〔公開〕

技術ノート 1（Eq. 57–58）：

$$
H^{\mathrm{grid},\sigma}\_{I\alpha,J\beta}(\mathbf L) = \Delta V \sum\_{\mathbf n}
\phi\_{I\alpha}(\mathbf r\_{\mathbf n} - \mathbf R\_I)\\, V^{\mathrm{grid}}\_\sigma(\mathbf r\_{\mathbf n})\\, \phi\_{J\beta}(\mathbf r\_{\mathbf n} - \mathbf R\_J - \mathbf L)
$$

（\\( \mathbf r\_{\mathbf n} \\) は周期的に拡張した格子点）。非共線では \\( \sigma \\) を 2×2 の成分 \\( (s,s') \\) に置き換えます。
この行列要素は、[疎ブロック形式](../design/sparse-block-format.md) の \\( H(\mathbf L) \\) の他の項（運動エネルギー、非局所、射影子展開した \\( V\_{\mathrm{na}} \\)）に加えます。

〔規約〕

- 和は 2.1 節の「2 つのカットオフ球の共通部分の格子点」だけで取ります。原子対は \\( \lVert\mathbf R\_J+\mathbf L-\mathbf R\_I\rVert \le r\_{c,I}+r\_{c,J} \\) のものです。
- 格子点上の基底関数の値は、動径関数表の補間で求め、単精度浮動小数点数で保持します。行列要素の相対誤差は \\( 10^{-7} \\) 程度になります。
  〔設計〕再実装は倍精度を既定とし、fixture との比較の許容誤差にこの差を見込みます。
- 粗い格子と細かい格子を使い分ける二重格子などの工夫はありません。格子は 1 種類で、密度・Poisson・XC・行列要素のすべてに使います
  （例外は 4.7 節の原子中心の球面積分）。

## 6. 力・応力への寄与（式のみ）〔公開〕〔規約〕

〔公開〕技術ノート 1（Eq. 63–69, 79–80）。格子に関係する力の項は、原子 \\( k \\) の位置 \\( \boldsymbol\tau\_k \\) について

$$
\frac{\partial E}{\partial\boldsymbol\tau\_k}\Big|\_{\mathrm{grid}}
{}= \sum\_\sigma\sum\_{\alpha\in k}\sum\_{J\beta,\mathbf L} 2\rho^\sigma\_{k\alpha,J\beta}(\mathbf L)\\,\Delta V\sum\_{\mathbf n}
\frac{\partial\phi\_{k\alpha}(\mathbf r\_{\mathbf n}-\boldsymbol\tau\_k)}{\partial\boldsymbol\tau\_k}\\, V^{\mathrm{grid}}\_\sigma(\mathbf r\_{\mathbf n})\\, \phi\_{J\beta}(\mathbf r\_{\mathbf n}-\mathbf R\_J-\mathbf L)
{}- \Delta V\sum\_{\mathbf n}\delta V\_H\\,\frac{\partial n^{(a)}\_k}{\partial\boldsymbol\tau\_k}
{}+ \Delta V\sum\_{\mathbf n}\sum\_\sigma v\_{\mathrm{xc},\sigma}\\,\frac{\partial n^{\mathrm{pcc}}\_{k,\sigma}}{\partial\boldsymbol\tau\_k}
$$

です（\\( n^{\mathrm{pcc}}\_{k,\sigma} \\) は 2.3 節でスピン \\( \sigma \\) に割り当てた分。密度行列の対称性で 2 倍）。
〔規約〕非共線計算の PCC の力は、局所スピン座標系でなく全体座標系の \\( v\_{\mathrm{xc}} \\) の成分を使うので、開殻内殻の種では厳密な微分と異なります（[非共線](noncollinear-energy-force.md) 3.4 節）。
基底関数の微分の項には \\( \rho \\) の変化（エネルギー密度行列の項）が別に付きます。これは [全エネルギーと力](energy-force.md) で扱います。

〔規約〕\\( E\_{\mathrm{xc}}^{(0)} \\) 補正（4.7 節）が on のとき、さらに

$$
{}-\Delta V\sum\_{\mathbf n} v^{\mathrm{CA}}\_{xc}(\hat n)\\,\frac{\partial\hat n\_k}{\partial\boldsymbol\tau\_k} + \frac{\partial E^{(0)}\_{\mathrm{xc}}}{\partial\boldsymbol\tau\_k},\qquad
\frac{\partial E^{(0)}\_{\mathrm{xc}}}{\partial\boldsymbol\tau\_k}
{}= \sum\_{I}\int \hat n\_I(\mathbf r-\mathbf R\_I)\\,\frac{d\varepsilon^{\mathrm{CA}}\_{xc}}{dn}(\hat n)\left[\frac{\partial\hat n\_k(\mathbf r-\boldsymbol\tau\_k)}{\partial\boldsymbol\tau\_k}(1-\delta\_{Ik}) - \delta\_{Ik}\sum\_{J\neq k}\frac{\partial\hat n\_J(\mathbf r-\mathbf R\_J)}{\partial\mathbf R\_J}\right] d^3r
$$

を加えます（\\( v^{\mathrm{CA}} = d(n\varepsilon^{\mathrm{CA}})/dn \\)。第 2 式の \\( I = k \\) の項は、積分領域が原子 \\( k \\) と一緒に動くことを
変数変換 \\( \mathbf u = \mathbf r - \boldsymbol\tau\_k \\) で取り込んだもの）。

応力：\\( E\_{\delta\mathrm H} \\) は 3.2 節の \\( \mathbf G \\) 空間の式を、\\( E\_{\mathrm{xc}} \\) は格子点の移動と \\( \Delta V \\) の変化、GGA の勾配項を通じて
格子ひずみで微分します。詳細は技術ノート 23 と [全エネルギーと力](energy-force.md) に譲ります（〔未確定〕）。

## 7. テスト

`tools/check_grid_poisson_xc.py` に、以下のうち解析的なものの参照実装があります。

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 分割数の規則：1.2 節の表の 4 例と、1.5 節の孤立系の箱（例：広がり 12, 12, 14.3 bohr、150 Ry → 辺 13.8515, 13.8515, 16.1601 bohr、\\( N = 54, 54, 63 \\)） | 厳密（整数一致） |
| 解析 | 周期 Poisson：歪んだ単位胞で \\( \delta n = \cos(\mathbf G\_0\cdot\mathbf r) + c \\) → \\( \delta V\_H = 4\pi\cos(\mathbf G\_0\cdot\mathbf r)/G\_0^2 \\)（定数 \\( c \\) は消える） | \\( 10^{-12} \\) |
| 解析 | 球形 Coulomb カットオフ：一辺 24 bohr の立方体、96³ 点、単位電荷の Gauss 分布 \\( (\alpha/\pi)^{3/2}e^{-\alpha r^2} \\)（\\( \alpha = 1.5 \\)）の中心から \\( r < L/2 - 4 \\) で \\( \delta V\_H = \mathrm{erf}(\sqrt\alpha r)/r \\)、自己エネルギー \\( \sqrt{\alpha/(2\pi)} \\) | \\( 10^{-8} \\) |
| 解析 | 同じ Gauss 分布をカットオフなし（3.1 節、平均を引く）で解くと自己エネルギーが 0.4296（孤立値 0.4886）になること：周期像の影響の大きさの確認 | \\( 10^{-4} \\) |
| 解析 | 27 点勾配：歪んだ格子で全次数 2 以下の多項式の勾配が厳密、直交格子で中心差分と一致 | \\( 10^{-10} \\) |
| 解析 | XC：下の参照値（libxc 7.1.2、非分極、Hartree）と \\( \varepsilon\_{xc} \\)・\\( v\_{xc} \\) を比較 | 相対 \\( 10^{-6} \\) |
| 解析 | 非共線：\\( \underline V\_{\mathrm{xc}} \\) の 12 成分の式、\\( \mathrm{Tr}(\underline V\delta\underline n) = \sum v\_\sigma\delta n'\_\sigma \\) | \\( 10^{-12} \\) |
| 不変量 | 電荷：\\( \Delta V\sum\_{\mathbf n} n \\) が電子数に \\( E \\) を上げると収束すること。\\( \Delta V\sum\delta n \\) は系の電荷の符号を変えた値に近い | 200 Ry で \\( 10^{-4} \\) 程度（要決定） |
| 不変量 | 格子原点の平行移動：全原子と `scf.fixed.grid` を同じベクトルだけずらすと、エネルギーが不変（格子点が原子と同じだけ動くため）。既定の原点（1.4 節）でも原点が原子の重心に付いて動くので同じ | 丸め誤差 |
| 解析 | 奇数の \\( N\_i \\) での波数の写し方（3.1 節）：\\( N=16 \\) で対称な選び方と一致、\\( N=15 \\) で最高波数の成分だけが異なること | \\( 10^{-12} \\) |
| 不変量 | 原点を \\( \mathbf h\_i \\) の分数だけずらしたときのエネルギーの変動（エッグボックス）が \\( E \\) とともに減り、\\( E\_{\mathrm{xc}}^{(0)} \\) 補正 on で off より小さいこと | 傾向の確認 |
| 不変量 | 力：格子項を含む力と、`scf.fixed.grid` を固定した全エネルギーの有限差分の一致 | \\( 10^{-4} \\) Ha/bohr 程度 |
| 照合 | 既存実装の出力（分割数、実際のカットオフ、格子原点、\\( E\_{\delta\mathrm H} \\)、\\( E\_{\mathrm{xc}} \\)、格子の行列要素）と fixture で比較。例：H₂O 孤立系（LDA・GGA-PBE）、Si 結晶（GGA-PBE）、bcc Fe（LSDA-CA、スピン分極）、非共線 Fe 鎖 | エネルギー \\( 10^{-6} \\) Ha、行列要素 \\( 10^{-6} \\) Ha（単精度の差を見込む。要決定） |

XC の参照値（libxc 7.1.2 で計算。`LDA,PZ` は Slater 交換＋PZ81、`LDA,PW` は Slater 交換＋PW92）：

| 汎関数 | 入力 | \\( \varepsilon\_{xc} \\) | \\( v\_{xc} \\) |
|---|---|---|---|
| PZ | \\( n = 10^{-4} \\) | −0.0495735123 | −0.0644958518 |
| PZ | \\( n = 10^{-2} \\) | −0.1970983191 | −0.2564000609 |
| PZ | \\( n = 1 \\) | −0.8091965677 | −1.0635669021 |
| PZ | \\( n = 10 \\)（\\( r\_s < 1 \\) の枝） | −1.6819531872 | −2.2215516600 |
| PW | \\( n = 10^{-2} \\) | −0.1968153660 | −0.2560329456 |
| PW | \\( n = 1 \\) | −0.8097590800 | −1.0642022422 |
| PZ（分極） | \\( n\_\uparrow = 0.08,\ n\_\downarrow = 0.02 \\) | −0.4167263619 | ↑ −0.5749254057、↓ −0.4326289466 |
| PW（分極） | \\( n\_\uparrow = 0.08,\ n\_\downarrow = 0.02 \\) | −0.4175794798 | ↑ −0.5772617587、↓ −0.4297856767 |
| PBE | \\( n = 0.1,\ \lvert\nabla n\rvert = 0.1 \\) | −0.3969182816 | \\( \partial e/\partial n \\) = −0.5149085316、\\( \partial e/\partial(\lvert\nabla n\rvert^2) \\) = −0.0156917755 |
| PBE | \\( n = 1,\ \lvert\nabla n\rvert = 0.5 \\) | −0.8097745831 | −1.0641456355、−0.0001245773 |

PBE の \\( v \\) は libxc の規約（\\( \partial e/\partial n \\) と \\( \partial e/\partial\sigma \\)、\\( \sigma = \lvert\nabla n\rvert^2 \\)）です。
\\( n = 1 \\) の PZ と PW の差（\\( 6\times10^{-4} \\) Ha）は、汎関数の取り違えを検出するのに十分な大きさです。

## 〔未確定〕

- 格子に関わる応力の各項の正確な式（技術ノート 23 との対応）。
- 照合 fixture の許容誤差（特に単精度の基底関数値による行列要素の差）の具体値。
- 電荷を持つ孤立系で `scf.coulomb.cutoff` を使わない場合の補正（Makov–Payne 型の補正の有無）。現時点では「補正しない」と理解していますが、エネルギー出力での確認が必要です。
