# 全エネルギーと力

> 状態：試作（レビュー待ち）。式は技術ノート 1（Ver. 1.3、以下 TN1）と
> 応力の技術ノート（Shiihara–Ozaki、以下 TN23）の方法を、このノートの記法で導出し直したものです。
> 数値で確かめた部分は `tools/check_energy_force.py` にあります。

## 記号〔設計〕

[記法](../design/notation.md) と [疎ブロック形式](../design/sparse-block-format.md) の記号に加えて、次を使います。

| 記号 | 意味 |
|---|---|
| \\( \sigma \\) | スピン（共線の場合 \\( \uparrow, \downarrow \\)） |
| \\( \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\) | 実空間の密度行列（スピンごと） |
| \\( \varepsilon^{\sigma}\_{I\alpha,J\beta}(\mathbf L) \\) | エネルギー重み付き密度行列（下で定義） |
| \\( O\_{I\alpha,J\beta}(\mathbf L) \\) | 原子対 \\( (I, J, \mathbf L) \\) の行列要素。\\( S \\)（重なり）、\\( T \\)（運動）、\\( V^{\mathrm{na}} \\)、\\( V^{\mathrm{NL}} \\) など |
| \\( \boldsymbol\tau\_{IJ\mathbf L} = \mathbf R\_J + \mathbf L - \mathbf R\_I \\) | 原子間ベクトル |
| \\( Z\_I \\) | 原子 \\( I \\) の価電子数（擬ポテンシャルの芯電荷） |
| \\( n\_\sigma(\mathbf r) \\)、\\( n = \sum\_\sigma n\_\sigma \\) | 価電子密度 |
| \\( n^{a}\_I(r) \\)、\\( n^{a} = \sum\_{I,\mathbf L} n^{a}\_I \\) | 孤立原子の価電子密度（PAO ファイルの密度）と、その重ね合わせ |
| \\( \delta n = n - n^{a} \\) | 差の密度 |
| \\( V^{a}\_{H,I}(r) \\) | \\( n^{a}\_I \\) が作る Hartree ポテンシャル。遠方で \\( Z\_I/r \\) |
| \\( \delta V\_H \\) | \\( \delta n \\) が作る Hartree ポテンシャル |
| \\( V^{\mathrm{na}}\_I = V^{\mathrm{loc}}\_I + V^{a}\_{H,I} \\) | 中性原子ポテンシャル（[中性原子](neutral-atom.md)） |
| \\( n^{\mathrm{pcc}} \\) | 部分内殻補正の密度（原子ごとの和） |
| \\( \mathbf r\_p \\)、\\( \Delta V \\) | 実空間格子点と、格子点 1 つあたりの体積 |
| \\( \Omega \\) | 単位胞の体積 |
| \\( r\_{c,I} \\) | 原子 \\( I \\) の PAO のカットオフ半径 |

原子対の和 \\( \sum\_{I,J,\mathbf L} \\) は、単位胞内の \\( I \\) と、\\( \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \le r\_{c,I} + r\_{c,J} \\) を満たす
すべての \\( (J, \mathbf L) \\)（\\( J = I \\)、\\( \mathbf L = 0 \\) を含む）にわたります。
以下では、和 \\( \sum\_{I,J,\mathbf L}\sum\_{\alpha\beta} \\) を \\( \sum\_{\mathrm{pairs}} \\) と略すことがあります。

### 密度行列とエネルギー重み付き密度行列〔公開〕

固有ベクトル \\( c^{\sigma}\_{n}(\mathbf k) \\)、固有値 \\( e^{\sigma}\_{n}(\mathbf k) \\)、占有数 \\( f^{\sigma}\_{n}(\mathbf k) \\)、
\\( \mathbf k \\) 点の重み \\( w\_{\mathbf k} \\)（和が 1）から

$$
\rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_n f^{\sigma}\_n(\mathbf k)\\, \mathrm{Re}\left[ e^{i\theta(\mathbf k,\mathbf L)}\\, c^{\sigma\ast}\_{I\alpha,n}(\mathbf k)\\, c^{\sigma}\_{J\beta,n}(\mathbf k) \right]
$$

$$
\varepsilon^{\sigma}\_{I\alpha,J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_n f^{\sigma}\_n(\mathbf k)\\, e^{\sigma}\_n(\mathbf k)\\, \mathrm{Re}\left[ e^{i\theta(\mathbf k,\mathbf L)}\\, c^{\sigma\ast}\_{I\alpha,n}(\mathbf k)\\, c^{\sigma}\_{J\beta,n}(\mathbf k) \right]
$$

とします（TN1 Eq. 4、74。位相 \\( \theta \\) は [疎ブロック形式](../design/sparse-block-format.md) の Bloch 和の規約に従います）。
共線スピンでない計算では、\\( \rho^{\uparrow} = \rho^{\downarrow} \\) をスピン 1 つ分として持ち、エネルギーでは 2 倍します〔規約〕。
以下の式は、共線（スピン分極なしを含む）で実数の基底・実数の \\( \rho \\) を前提にします。非共線は〔未確定〕です。

## 全エネルギーの分割〔公開〕

TN1 Eq. 7 の素朴な分割は

$$
E = E\_{\mathrm{kin}} + E^{\mathrm{loc}}\_{\mathrm{ec}} + E^{\mathrm{NL}}\_{\mathrm{ec}} + E\_{\mathrm{ee}} + E\_{\mathrm{xc}} + E\_{\mathrm{cc}}
$$

です。\\( E^{\mathrm{loc}}\_{\mathrm{ec}} = \int n \sum\_I V^{\mathrm{loc}}\_I \\)、
\\( E\_{\mathrm{ee}} = \frac12\iint n(\mathbf r)n(\mathbf r')/\lvert\mathbf r-\mathbf r'\rvert \\)、
\\( E\_{\mathrm{cc}} = \frac12\sum\_{I\ne J} Z\_I Z\_J/\lvert\mathbf R\_I-\mathbf R\_J\rvert \\)（周期像を含む）です。
長距離の項どうしの打ち消しを明示するため、\\( n = n^{a} + \delta n \\) を代入して並べ替えます。

$$
E\_{\mathrm{ee}} = \tfrac12\int n^{a} V^{a}\_H + \int \delta n\\, V^{a}\_H + \tfrac12\int \delta n\\, \delta V\_H,\qquad
\int n\\, V^{\mathrm{na}} = E^{\mathrm{loc}}\_{\mathrm{ec}} + \int n^{a} V^{a}\_H + \int \delta n\\, V^{a}\_H
$$

から（\\( V^{\mathrm{na}} = \sum\_I V^{\mathrm{na}}\_I \\)、\\( V^{a}\_H = \sum\_I V^{a}\_{H,I} \\)）、

$$
E^{\mathrm{loc}}\_{\mathrm{ec}} + E\_{\mathrm{ee}} + E\_{\mathrm{cc}}
{}= \underbrace{\int n\\, V^{\mathrm{na}}}\_{E\_{\mathrm{na}}}
{}+ \underbrace{\tfrac12\int \delta n\\, \delta V\_H}\_{E\_{\delta\mathrm{ee}}}
{}+ \underbrace{E\_{\mathrm{cc}} - \tfrac12\int n^{a} V^{a}\_H}\_{E\_{\mathrm{scc}}}
$$

が得られます（TN1 Eq. 19–22）。したがって

$$
E = E\_{\mathrm{kin}} + E\_{\mathrm{na}} + E^{\mathrm{NL}}\_{\mathrm{ec}} + E\_{\delta\mathrm{ee}} + E\_{\mathrm{xc}} + E\_{\mathrm{scc}}
$$

です（TN1 Eq. 23、TN23 Eq. 1）。

- \\( E\_{\mathrm{scc}} = \frac12\sum\_{I,J}\left[ \frac{Z\_I Z\_J}{\lvert\cdot\rvert}(1-\delta\_{IJ}^{\mathbf L=0}) - \int n^{a}\_I V^{a}\_{H,J} \right] \\) の各対は、
  \\( n^{a}\_I \\) と \\( n^{a}\_J \\) の台が重ならない距離では厳密に 0 になります
  （球対称な電荷 \\( Z\_I \\) の外側のポテンシャルは点電荷と同じなので）。つまり遮蔽された短距離の対ポテンシャルです。
  和の \\( I = J \\)、\\( \mathbf L = 0 \\) の項は原子の自己エネルギー \\( -\frac12\int n^{a}\_I V^{a}\_{H,I} \\) で、配置によらない定数です。
- \\( E\_{\mathrm{na}} \\) の \\( V^{\mathrm{na}}\_I \\) は中性原子の外で 0 なので、行列要素は有限の原子対だけを含みます。
- 長距離の静電相互作用は \\( E\_{\delta\mathrm{ee}} \\) だけが担い、これは FFT で解きます（[グリッド](grid-poisson-xc.md)）。
- 交換相関は部分内殻補正を含めて
  \\( E\_{\mathrm{xc}} = \int (n + n^{\mathrm{pcc}})\\, \epsilon\_{\mathrm{xc}}(n\_\uparrow + \tfrac12 n^{\mathrm{pcc}}, n\_\downarrow + \tfrac12 n^{\mathrm{pcc}}) \\) です（TN1 Eq. 15。GGA では勾配も引数に入る）。

## 出力されるエネルギー成分〔規約〕

出力ファイル（`System.Name.out`）の「Total energy (Hartree)」の節には、次のラベルで値が出ます。
ラベルの後ろにピリオドが付きます（`Utot.` など）。単位は Hartree です。

| ラベル | 定義（このノートの記法） | 評価の方法 |
|---|---|---|
| `Uele` | \\( \sum\_\sigma\sum\_{\mathbf k} w\_{\mathbf k}\sum\_n f^{\sigma}\_n e^{\sigma}\_n \\)（バンドエネルギー。`Utot` には入らない） | 固有値の和 |
| `Ukin` | \\( \sum\_\sigma \sum\_{\mathrm{pairs}} \rho^{\sigma} T \\) | 行列のトレース（二中心積分） |
| `UH0` | \\( -\frac12\sum\_{I,J,\mathbf L} s\_{IJ} \int n^{a}\_I V^{a}\_{H,J} \\)（自己項 \\( J=I, \mathbf L=0 \\) を含む。\\( s\_{IJ} \\) は下記） | 原子対ごとの数値積分 |
| `UH1` | \\( E\_{\delta\mathrm{ee}} = \frac12\Delta V\sum\_p \delta n(\mathbf r\_p)\\, \delta V\_H(\mathbf r\_p) \\) | 格子上の和 |
| `Una` | \\( E\_{\mathrm{na}} \\) | 行列のトレース（射影子展開）または格子上の和 |
| `Unl` | \\( \sum\_\sigma \sum\_{\mathrm{pairs}} \rho^{\sigma} V^{\mathrm{NL}} \\) | 行列のトレース（二中心積分の積） |
| `Uxc0`、`Uxc1` | \\( E\_{\mathrm{xc}} \\) のスピンごとの部分（下記） | 格子上の和 ＋ 原子密度の補正 |
| `Ucore` | \\( \frac12\sum\_{I,(J,\mathbf L)\ne(I,0)} Z\_I Z\_J / \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \\) | 原子対の和 |
| `Uhub` | DFT+U のエネルギー | 〔未確定〕 |
| `Ucs` | スピンの向きの拘束エネルギー | 〔未確定〕 |
| `Uzs`、`Uzo` | スピン・軌道磁気モーメントの Zeeman 項 | 〔未確定〕 |
| `Uef` | 外部電場のエネルギー \\( \Delta V\sum\_p n(\mathbf r\_p) V\_{\mathrm{ef}}(\mathbf r\_p) \\)（電子のみ） | 格子上の和 |
| `UvdW` | DFT-D2 または D3 の分散エネルギー | 原子対（三体）の和 |
| `Uch` | 内殻正孔を作るためのペナルティ項 | 〔未確定〕 |
| `Utot` | 下の和 | |
| `UpV` | 加えた圧力と体積の積（`MD.applied.pressure`） | |
| `Enpy` | `Utot` ＋ `UpV`（エンタルピー） | |

全エネルギーは

$$
U\_{\mathrm{tot}} = U\_{\mathrm{kin}} + U\_{\mathrm{H0}} + U\_{\mathrm{H1}} + U\_{\mathrm{na}} + U\_{\mathrm{nl}} + U\_{\mathrm{xc0}} + U\_{\mathrm{xc1}} + U\_{\mathrm{core}} + U\_{\mathrm{hub}} + U\_{\mathrm{cs}} + U\_{\mathrm{zs}} + U\_{\mathrm{zo}} + U\_{\mathrm{ef}} + U\_{\mathrm{vdW}} + U\_{\mathrm{ch}}
$$

です。TN1 の記号との対応は \\( E\_{\mathrm{scc}} \leftrightarrow U\_{\mathrm{core}} + U\_{\mathrm{H0}} \\)、
\\( E\_{\delta\mathrm{ee}} \leftrightarrow U\_{\mathrm{H1}} \\)、\\( E\_{\mathrm{xc}} \leftrightarrow U\_{\mathrm{xc0}} + U\_{\mathrm{xc1}} \\) です。

- 電子温度があっても、`Utot` には電子のエントロピー項 \\( -TS \\) を含めません。
- スピン分極なしでは `Uxc0` と `Uxc1` は等しく、それぞれ \\( E\_{\mathrm{xc}} \\) の半分です。
  共線スピンでは \\( U\_{\mathrm{xc}\sigma} = \Delta V\sum\_p (n\_\sigma + \tfrac12 n^{\mathrm{pcc}})\\, \epsilon\_{\mathrm{xc}}(\mathbf r\_p) \\)
  （\\( \epsilon\_{\mathrm{xc}} \\) は全密度で決まる 1 電子あたりのエネルギー）に、下の補正の半分ずつを加えたものです。
- `scf.EigenvalueSolver` が `Cluster-LNO` のときは、別の式（後述の Harris 型）で `Utot` を作り、
  出力のラベルも `UH_dc`、`Uxc0_dc`、`Uxc1_dc` を含む別の並びになります。

### 行列のトレースで評価する項〔規約〕

`Ukin`、`Una`（射影子展開のとき）、`Unl` は、密度行列と行列要素の、保持している原子対の範囲での積和です。

$$
U\_X = g\_s \sum\_{\sigma} \sum\_{I,J,\mathbf L}\sum\_{\alpha\beta} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, X\_{I\alpha,J\beta}(\mathbf L)
$$

\\( g\_s \\) はスピン分極なしで 2、それ以外で 1 です。
\\( T \\) は [二中心積分](two-center-overlap-kinetic.md)、\\( V^{\mathrm{NL}} \\) は [非局所擬ポテンシャル](nonlocal-pseudopotential.md)、
\\( V^{\mathrm{na}} \\) は [中性原子](neutral-atom.md) の章で作ります。
`scf.ProExpn.VNA`（既定 on）が off のときは、`Una` を格子上の和 \\( \Delta V\sum\_p n(\mathbf r\_p) V^{\mathrm{na}}(\mathbf r\_p) \\) で評価します。

### 格子上で評価する項〔規約〕

格子点の値 \\( n\_\sigma(\mathbf r\_p) \\)、\\( n^{a}(\mathbf r\_p) \\)、\\( n^{\mathrm{pcc}}(\mathbf r\_p) \\)、\\( \delta V\_H(\mathbf r\_p) \\)、
\\( \epsilon\_{\mathrm{xc}}(\mathbf r\_p) \\) から、台形則（格子点の単純和に \\( \Delta V \\) を掛ける）で積分します。

- `UH1`：\\( \frac12\Delta V\sum\_p (n - n^{a})\\, \delta V\_H \\)。
- `Uxc`：\\( \Delta V\sum\_p \sum\_\sigma (n\_\sigma + \frac12 n^{\mathrm{pcc}})\\, \epsilon\_{\mathrm{xc}} \\)。
  格子点上の \\( \epsilon\_{\mathrm{xc}} \\) は、`Uele` の Hamiltonian と同じ汎関数で、同じ格子・同じ勾配の差分（GGA の場合）で求めます。

### 交換相関の原子密度補正〔規約〕

`scf.Exc0.correction.flag`（既定 on）が on のとき、`Uxc0 + Uxc1` に次の \\( \Delta E^{(0)} \\) を加えます（半分ずつ）。

$$
\Delta E^{(0)} = -\Delta V\sum\_p \left(n^{a} + n^{\mathrm{pcc}}\right)\epsilon^{\mathrm{CA}}\\!\left(n^{a} + n^{\mathrm{pcc}}\right)\Big\rvert\_{\mathbf r\_p}
{}+ \sum\_{I}\int\_{\lVert\mathbf x\rVert < r\_{c,I}} d^3x\\; n^{a}\_I(\mathbf x)\\, \epsilon^{\mathrm{CA}}\\!\left(n^{a}(\mathbf R\_I + \mathbf x)\right)
$$

- \\( \epsilon^{\mathrm{CA}} \\) は、計算で選んだ汎関数にかかわらず、常にスピン分極なしの Ceperley–Alder LDA
  です（パラメータ化の詳細は [グリッド](grid-poisson-xc.md) の章で確定させる）。
- 第 1 項は粗い実空間格子上の和、第 2 項は原子中心の細かい求積です。原子 \\( I \\) を中心とする球の中で、
  動径方向は \\( [0, r\_{c,I}] \\) の Gauss–Legendre 150 点、角度方向は Lebedev 則
  （点数は `scf.Num.Leb.Grid`、既定 590）を使います。
  \\( n^{a}(\mathbf R\_I + \mathbf x) \\) は、近接原子（原子対の範囲にある原子）の \\( n^{a}\_J \\) の和で、各 \\( n^{a}\_J \\) は自分のカットオフの外で 0 とします。
- 意図は、原子密度の重ね合わせの交換相関エネルギーを、格子でなく精度のよい原子中心の求積で評価し直すことです
  （格子の粗さによる誤差、とくに格子に対する原子の位置で変わる誤差を減らす）。
- 第 2 項の密度には \\( n^{\mathrm{pcc}} \\) が含まれず、第 1 項には含まれます。部分内殻補正のある原子では
  2 つの項は同じ量の 2 通りの評価になっておらず、差 \\( \int[(n^{a}+n^{\mathrm{pcc}})\epsilon(n^{a}+n^{\mathrm{pcc}}) - n^{a}\epsilon(n^{a})] \\) が残ります。
  既存実装の挙動として記録します。再実装でもこの定義に合わせるかどうかは〔未確定〕です。

### 遮蔽された原子間項：Ucore と UH0〔公開〕〔規約〕

〔公開〕TN1 は \\( E\_{\mathrm{scc}} \\) の第 2 項が原子間距離だけの関数なので、数値的に表にして補間すると述べています（TN1 Eq. 22 の後）。

〔規約〕既存実装は、距離の表ではなく、原子対ごとに次の数値積分を毎回評価します。

$$
J\_{IJ}(d) = \int d^3 r\\; n^{a}\_I(\mathbf r)\\, V^{a}\_{H,J}(\mathbf r - d\\,\hat{\mathbf z})
$$

- 積分は対の軸を \\( z \\) とする円筒座標で行います。\\( n^{a}\_I \\) の台 \\( [-r\_{c,I}, r\_{c,I}] \\) を \\( z \\) 方向に等間隔に分け
  （間隔の目安 \\( \pi/\sqrt{9600} \approx 0.032 \\) bohr。両端を含む奇数個の点になるよう間隔を調整）、
  円筒の半径方向 \\( [0, r\_{c,I}] \\) は Gauss–Legendre 150 点、方位角は解析的に \\( 2\pi \\) です。
  \\( n^{a}\_I \\) は PAO ファイルの原子密度、\\( V^{a}\_{H,J} \\) は擬ポテンシャルファイルの原子 Hartree ポテンシャル
  （範囲外では \\( Z\_J/r \\)）を動径方向に補間します。
- \\( U\_{\mathrm{H0}} \\) は、各対について \\( J\_{IJ} \\) と \\( J\_{JI} \\) の平均をとり、\\( -\frac12 \\) を掛けて和をとります（対称化）。
- 原子の異なる対（\\( J \ne I \\) または \\( \mathbf L \ne 0 \\)）には、種の組ごとの定数

$$
s\_{IJ} = \frac{Z\_I Z\_J / (r\_{c,I} + r\_{c,J})}{J\_{IJ}(r\_{c,I} + r\_{c,J})}
$$

  を掛けます（自己項は 1）。これにより、原子対の範囲の端 \\( d = r\_{c,I} + r\_{c,J} \\) で
  \\( Z\_I Z\_J/d - s\_{IJ} J\_{IJ}(d) = 0 \\) が数値的にも成り立ち、対が範囲に出入りしてもエネルギーが連続になります。
  \\( s\_{IJ} \\) は最初の構造で一度だけ求めて固定し、力では微分しません。
- \\( U\_{\mathrm{core}} \\) は同じ原子対の範囲での点電荷の和です。範囲の外の対は \\( U\_{\mathrm{core}} \\) にも \\( U\_{\mathrm{H0}} \\) にも入りませんが、
  遮蔽された対の和は範囲の外で 0 なので、合計は打ち切りの影響を受けません。
  一方で、\\( U\_{\mathrm{core}} \\) と \\( U\_{\mathrm{H0}} \\) を別々に見ると、どちらも打ち切りに依存する大きな数です。比較は和で行います。

〔設計〕再実装では、\\( J\_{IJ}(d) \\) を種の組ごとに \\( d \\) の表として前計算し、スプラインで補間しても構いません
（TN1 の方法）。その場合も、端点での規格化 \\( s\_{IJ} \\) と、表の補間と整合する微分を使います。
円筒格子の求積は、Gauss 型の密度と、その Hartree ポテンシャル \\( Z\\,\mathrm{erf}(\sqrt{a}r)/r \\) の組で、
解析解 \\( Z\_I Z\_J\\,\mathrm{erf}(\sqrt{\mu}d)/d \\)（\\( \mu = a\_I a\_J/(a\_I+a\_J) \\)）と相対 \\( 10^{-10} \\) で一致することを確認しています。

### その他の項

- 〔公開〕`UvdW`：Grimme の DFT-D2・D3（マニュアル 54 節）。原子位置だけの関数なので、力・応力は解析微分で足します。
- 〔未確定〕`Uhub`（DFT+U の占有数の定義：on-site、full、dual）、`Ucs`、`Uzs`、`Uzo`、`Uch` の定義と力。
  これらは対応する機能の章で扱います。
- 〔規約〕`Uef` は電子密度と電場ポテンシャルの積の格子和で、原子核（芯電荷）と電場の相互作用は含みません。

## 整合性の恒等式

### バンドエネルギーと二重計数〔規約〕

Hamiltonian の行列要素は

$$
H^{\sigma} = T + V^{\mathrm{na}} + V^{\mathrm{NL}} + \langle \phi \rvert \delta V\_H + v\_{\mathrm{xc},\sigma} \lvert \phi \rangle\_{\mathrm{grid}} + \cdots
$$

（最後の項は格子上の和で作る行列要素）なので、同じ固有ベクトルから作った \\( \rho \\) について
\\( U\_{\mathrm{ele}} = g\_s\sum\_\sigma\sum\_{\mathrm{pairs}} \rho^{\sigma} H^{\sigma} \\) が恒等的に成り立ちます。
これを使うと、SCF が収束したとき

$$
U\_{\mathrm{tot}} = U\_{\mathrm{ele}}
{}- \Delta V\sum\_p n\\,\delta V\_H + U\_{\mathrm{H1}}
{}- \Delta V\sum\_p \sum\_\sigma n\_\sigma v\_{\mathrm{xc},\sigma} + U\_{\mathrm{xc0}} + U\_{\mathrm{xc1}}
{}+ U\_{\mathrm{core}} + U\_{\mathrm{H0}} + \cdots
$$

（Harris 型、二重計数の形）が成り立ちます。`Cluster-LNO` ではこの形で `Utot` を作り、
\\( -\Delta V\sum n\\,\delta V\_H \\) を `UH_dc`、\\( -\Delta V\sum n\_\sigma v\_{\mathrm{xc},\sigma} \\) を `Uxc0_dc`、`Uxc1_dc` として出力します。
二重計数の項は価電子密度 \\( n\_\sigma \\) だけを使います（\\( n^{\mathrm{pcc}} \\) は Hamiltonian の固有値に寄与しないため）。
二重計数の Hartree 項が \\( \delta n \\) でなく \\( n \\) を使うのは、Hamiltonian の \\( \delta V\_H \\) が全価電子に作用するからです。

SCF が途中のとき（Hamiltonian を作った密度と、その固有ベクトルから作った密度が異なるとき）は 2 つの式は一致しません。
テストでは十分に収束させてから比べます。

### その他の恒等式〔設計〕

- 孤立原子 1 個（周期像が原子対の範囲に入らない大きな箱）：\\( U\_{\mathrm{core}} = 0 \\)、
  \\( U\_{\mathrm{H0}} = -\frac12\int n^{a}\_I V^{a}\_{H,I} \\)。
- 原子対の範囲の端：遮蔽された対のエネルギーが \\( s\_{IJ} \\) の定義により 0。
- \\( \Delta V\sum\_p \delta n(\mathbf r\_p) \approx 0 \\)（系が中性なので。格子の誤差の範囲で）。

## 力〔公開〕

### 符号と単位〔規約〕

- 力は \\( \mathbf F\_K = -\partial U\_{\mathrm{tot}}/\partial\mathbf R\_K \\)（Hartree/bohr）です。
  出力ファイルの `<coordinates.forces` ブロックには、原子番号、元素記号、座標（Å）、力の 3 成分（Hartree/bohr）が並びます。
  MD の履歴ファイル（`.md`）も同じ符号の力を出します。
- 構造最適化と `MD.Type NOMD` では、力に総和を 0 にする補正をかけません。
  分子動力学（`NVE`、`NVT_VS`、`NVT_VS2`、`NVT_VS4`、`NVT_NH`、`NVT_Langevin`）では、時間発展の前に、
  固定されていない原子について成分ごとに力の平均を引き、総和を 0 にします。
  その場合に出力ファイルの力が補正前か補正後かは〔未確定〕です。
- 力は SCF の最後の密度行列で、全エネルギーの解析微分として評価します。

### 全体の形〔公開〕

\\( \mathbf R\_K \\) に関する微分は、TN1 Eq. 59–81 にならって次のように分けられます。
密度行列の変化を通じた寄与は、固有ベクトルの正規直交条件 \\( c^\dagger S c = 1 \\) の微分
\\( \partial c^\dagger S c + c^\dagger S \partial c = -c^\dagger (\partial S) c \\)（TN1 Eq. 73）を使うと、すべて重なり行列の微分にまとまります
（Pulay の項）。

$$
\frac{\partial U\_{\mathrm{tot}}}{\partial \mathbf R\_K}
{}= \underbrace{- g\_s\sum\_\sigma\sum\_{\mathrm{pairs}} \varepsilon^{\sigma}\\, \frac{\partial S}{\partial\mathbf R\_K}}\_{\text{Pulay}}
{}+ \underbrace{g\_s\sum\_\sigma\sum\_{\mathrm{pairs}} \rho^{\sigma} \frac{\partial (T + V^{\mathrm{na}} + V^{\mathrm{NL}})}{\partial\mathbf R\_K}}\_{\text{行列要素}}
{}+ \underbrace{\frac{\partial U\_{\mathrm{grid}}}{\partial\mathbf R\_K}\Big\rvert\_{\rho}}\_{\text{格子}}
{}+ \frac{\partial (U\_{\mathrm{core}} + U\_{\mathrm{H0}})}{\partial\mathbf R\_K}
{}+ \cdots
$$

最初の 2 項の和は、\\( \rho \\) を固定した行列要素の微分と、Pulay の項で、
\\( \sum f e \\) の微分の公式 \\( \partial \sum\_n f\_n e\_n = \mathrm{Tr}[\rho\\, \partial H] - \mathrm{Tr}[\varepsilon\\, \partial S] \\)（\\( f \\) 固定）に当たります。
この公式は `tools/check_energy_force.py` で一般化固有値問題について数値的に確認しています。

### 行列要素の微分〔公開〕

- 二中心積分 \\( O(\boldsymbol\tau) \\)（\\( S \\)、\\( T \\)）は \\( \boldsymbol\tau = \mathbf R\_J + \mathbf L - \mathbf R\_I \\) だけの関数なので、
  \\( \partial O/\partial\mathbf R\_I = -\nabla\_{\boldsymbol\tau} O \\)、\\( \partial O/\partial\mathbf R\_J = +\nabla\_{\boldsymbol\tau} O \\)、それ以外の原子については 0 です。
  \\( \nabla\_{\boldsymbol\tau} O \\) は [二中心積分](two-center-overlap-kinetic.md) の 6 節の解析微分を使います。
  同じ原子上の対（\\( \boldsymbol\tau = 0 \\) で \\( I = J \\)）は配置によらないので寄与しません。
- 非局所項 \\( V^{\mathrm{NL}}\_{I\alpha,J\beta} = \sum\_K\sum\_{ij}\langle\phi\_{I\alpha}\vert\beta^K\_i\rangle D^K\_{ij}\langle\beta^K\_j\vert\phi\_{J\beta}\rangle \\) は
  3 つの原子 \\( I, J, K \\) の位置に依存します。各因子は二中心積分なので、

$$
\frac{\partial V^{\mathrm{NL}}\_{I\alpha,J\beta}}{\partial \mathbf R\_K}
{}= \sum\_{ij}\left[\frac{\partial\langle\phi\_{I\alpha}\vert\beta^K\_i\rangle}{\partial\mathbf R\_K} D^K\_{ij}\langle\beta^K\_j\vert\phi\_{J\beta}\rangle
{}+ \langle\phi\_{I\alpha}\vert\beta^K\_i\rangle D^K\_{ij}\frac{\partial\langle\beta^K\_j\vert\phi\_{J\beta}\rangle}{\partial\mathbf R\_K}\right]
$$

  のように、因子ごとに二中心積分の微分を使います。\\( K \\) が \\( I \\) や \\( J \\) と一致する場合も、因子ごとに同じ規則を当てはめます。
  並進不変性から \\( \partial/\partial\mathbf R\_I + \partial/\partial\mathbf R\_J + \partial/\partial\mathbf R\_K = 0 \\) が成り立つので、
  2 つを計算すれば 3 つ目が得られます（TN1 Eq. 62）。
- 中性原子項を射影子展開で作る場合（TN1 Eq. 32–34、[中性原子](neutral-atom.md)）も、
  \\( \langle\phi\vert V^{\mathrm{na}}\_K \bar R\_{l\zeta} Y\_{lm}\rangle \\) を二中心積分として同じ形で微分します。
  展開を使わずに直接計算する組み合わせ（一中心・二中心の場合）の微分は〔未確定〕です。

### 格子上の項の微分〔公開〕

格子点 \\( \mathbf r\_p \\) は原子と一緒に動かず固定されたままで、原子に付随する関数（基底関数、\\( n^{a}\_K \\)、\\( n^{\mathrm{pcc}}\_K \\)、
格子モードの \\( V^{\mathrm{na}}\_K \\)）だけが動くとして、離散化したエネルギーを厳密に微分します（TN1 Eq. 63–69、79–80）。
このため、解析的な力は任意の格子間隔で、離散化した全エネルギーの有限差分と一致します。

関数 \\( g(\mathbf r\_p - \mathbf R\_K) \\) の微分は \\( \partial g/\partial\mathbf R\_K = -\nabla g(\mathbf r\_p - \mathbf R\_K) \\) です。

1. **基底関数の移動**（TN1 Eq. 67、79）：\\( \rho \\) を固定したときの密度の変化

$$
\frac{\partial n\_\sigma(\mathbf r\_p)}{\partial\mathbf R\_K}\Big\rvert\_{\rho}
{}= -2 \sum\_{\alpha}\sum\_{J,\mathbf L,\beta} \rho^{\sigma}\_{K\alpha,J\beta}(\mathbf L)\\, \nabla\phi\_{K\alpha}(\mathbf r\_p - \mathbf R\_K)\\, \phi\_{J\beta}(\mathbf r\_p - \mathbf R\_J - \mathbf L)
$$

   （係数 2 は \\( \rho \\) の対称性から \\( K \\) が左右どちらに現れる場合もまとめたもの）に、
   有効ポテンシャル \\( V^{\sigma}\_{\mathrm{eff}} = \delta V\_H + v\_{\mathrm{xc},\sigma} \\)（`Una` を格子で評価するときは \\( + V^{\mathrm{na}} \\)、
   電場があるときは \\( + V\_{\mathrm{ef}} \\)）を掛けて \\( g\_s \Delta V\sum\_p\sum\_\sigma \\) で和をとります。
   \\( \delta V\_H \\) を通じた間接の変化は、\\( E\_{\delta\mathrm{ee}} \\) が \\( \delta n \\) の 2 次形式であることから、この 1 項にまとまります（TN1 Eq. 64–65）。
2. **原子密度の移動**（TN1 Eq. 66、68、80）：\\( E\_{\delta\mathrm{ee}} \\) は \\( n^{a} \\) にも \\( \delta n = n - n^{a} \\) を通じて依存するので、
   \\( -\Delta V\sum\_p \delta V\_H(\mathbf r\_p)\\, \partial n^{a}\_K(\mathbf r\_p - \mathbf R\_K)/\partial\mathbf R\_K \\)。
3. **部分内殻補正の移動**（TN1 Eq. 69、80）：\\( \frac12\Delta V\sum\_p\sum\_\sigma v\_{\mathrm{xc},\sigma}(\mathbf r\_p)\\, \partial n^{\mathrm{pcc}}\_K/\partial\mathbf R\_K \\)。
4. **格子モードの中性原子ポテンシャルの移動**：\\( \Delta V\sum\_p n(\mathbf r\_p)\\, \partial V^{\mathrm{na}}\_K(\mathbf r\_p - \mathbf R\_K)/\partial\mathbf R\_K \\)。
5. **原子密度の補正 \\( \Delta E^{(0)} \\) の微分**〔規約〕：第 1 項は
   \\( -\Delta V\sum\_p v^{\mathrm{CA}}(n^{a}+n^{\mathrm{pcc}})\\, \partial(n^{a}\_K + n^{\mathrm{pcc}}\_K)/\partial\mathbf R\_K \\)
   （\\( v^{\mathrm{CA}} = d(n\epsilon^{\mathrm{CA}})/dn \\)）。第 2 項は原子中心の求積点が中心原子と一緒に動くので、
   中心 \\( I \\) と近接原子 \\( J \\) の相対位置だけに依存し、

$$
\frac{\partial}{\partial \mathbf R\_J}\int n^{a}\_I(\mathbf x)\\, \epsilon^{\mathrm{CA}}(n^{a}(\mathbf R\_I + \mathbf x))\\, d^3x
{}= -\int n^{a}\_I(\mathbf x)\\, \frac{d\epsilon^{\mathrm{CA}}}{dn}\\, \nabla n^{a}\_J(\mathbf R\_I + \mathbf x - \mathbf R\_J)\\, d^3x
$$

   を \\( J \\) に、その符号を反転したものを \\( I \\) に加えます（同じ求積点で評価）。

GGA の場合、\\( v\_{\mathrm{xc}} \\) は離散化したエネルギー \\( \Delta V\sum\_p f(n\_p, \nabla\_h n\_p) \\)（\\( \nabla\_h \\) は下の応力の節で述べる差分勾配）を
\\( n\_p \\) で偏微分したもの（差分勾配の転置を含む）でなければ、上の 1 の式は有限差分と一致しません。
既存実装の \\( v\_{\mathrm{xc}} \\) がこの意味で離散的に厳密かどうかは〔未確定〕です。

### 遮蔽された原子間項の微分〔規約〕

- \\( U\_{\mathrm{core}} \\)：\\( \partial/\partial\mathbf R\_I \\) は \\( -\sum\_{(J,\mathbf L)} Z\_I Z\_J\\, (\mathbf R\_I - \mathbf R\_J - \mathbf L)/d^3 \\)（\\( d = \lVert\boldsymbol\tau\rVert \\)）。
- \\( U\_{\mathrm{H0}} \\)：\\( J\_{IJ}(d) \\) の \\( d \\) による微分を、同じ円筒格子の上で \\( V^{a}\_{H,J} \\) の動径微分の補間から評価し、
  対の軸方向の単位ベクトルを掛けます。係数 \\( s\_{IJ} \\) は定数として扱います。
  \\( J\_{IJ} \\) の微分は \\( n^{a}\_I \\) を固定して \\( V^{a}\_{H,J} \\) を動かした形で求め、\\( J\_{JI} \\) についても同様にして、両者を対称化します。

### 電子温度と力〔設計〕

上の解析的な力は、占有数 \\( f \\) を Fermi–Dirac 分布にしても、エントロピー項を含む自由エネルギー
\\( A = U\_{\mathrm{tot}} - T\_e S \\)（\\( S = -g\_s\sum w\_{\mathbf k}\sum\_n [f\ln f + (1-f)\ln(1-f)] \\)）の微分になります
（\\( f \\) についての停留性と電子数の保存から）。`Utot` 自身の微分ではないので、
金属や高い電子温度では、`Utot` の有限差分と解析的な力が \\( O(T\_e) \\) で食い違います（TN1 の注意と同じ趣旨）。
力のテストでは、絶縁体を使うか、\\( A \\) の有限差分と比べます。

## 応力〔公開〕

### 定義〔公開〕

TN23 に従い、ひずみ \\( \boldsymbol\epsilon \\)（対称な 3×3 行列）のもとで、位置ベクトルを \\( \mathbf r \to (1+\boldsymbol\epsilon)\mathbf r \\)、
格子ベクトルと原子位置も同じ変換で動かし、基底関数・射影子・原子密度などの形は変えない（中心だけが動く）とします（TN23 Eq. 5–9）。
エネルギーの 1 次の係数

$$
\Sigma\_{\gamma\eta} = \frac{\partial U\_{\mathrm{tot}}}{\partial\epsilon\_{\gamma\eta}}\Big\rvert\_{\boldsymbol\epsilon = 0}
$$

を求めます（TN23 Eq. 10–11）。逆格子ベクトルは \\( (1+\boldsymbol\epsilon)^{-T} \\) で変わり、\\( \mathbf k\cdot\mathbf L \\)、\\( \mathbf G\cdot\mathbf r \\) は不変です。
単位胞の積分は \\( \det(1+\boldsymbol\epsilon) \\) 倍になるので、格子上の和には \\( \delta\_{\gamma\eta}\times \\)（その項）の「体積項」が出ます。

原子位置の関数 \\( g(\boldsymbol\tau) \\) のひずみ微分は、連鎖律で

$$
\frac{\partial g}{\partial\epsilon\_{\gamma\eta}} = \frac{\partial g}{\partial\tau\_\gamma}\\,\tau\_\eta
$$

です（TN23 Eq. 14–15）。\\( \boldsymbol\epsilon \\) は対称なので、最後に \\( \Sigma \\) を対称化します〔規約〕。

### 各項の寄与〔公開〕

| 項 | \\( \partial/\partial\epsilon\_{\gamma\eta} \\) | 出典 |
|---|---|---|
| 行列要素（\\( T \\)、\\( V^{\mathrm{NL}} \\)、射影子展開の \\( V^{\mathrm{na}} \\)） | \\( g\_s\sum\_\sigma\sum\_{\mathrm{pairs}}\rho^{\sigma}\sum\_{\text{因子}} \frac{\partial O}{\partial\tau\_\gamma}\tau\_\eta \\)。三中心の項は二中心の因子ごとに、その因子の原子間ベクトルを使う | TN23 Eq. 16–31 |
| Pulay（重なり） | \\( -g\_s\sum\_\sigma\sum\_{\mathrm{pairs}}\varepsilon^{\sigma}\frac{\partial S}{\partial\tau\_\gamma}\tau\_\eta \\) | TN23 Eq. 62–67 |
| \\( E\_{\delta\mathrm{ee}} \\) | \\( \delta\_{\gamma\eta}E\_{\delta\mathrm{ee}} + \int \frac{\partial\delta n}{\partial\epsilon\_{\gamma\eta}}\delta V\_H + 4\pi\Omega\sum\_{\mathbf G\ne 0}\lvert\delta n(\mathbf G)\rvert^2\frac{G\_\gamma G\_\eta}{G^4} \\) | TN23 Eq. 32–45 |
| \\( E\_{\mathrm{xc}} \\)（LDA） | \\( \delta\_{\gamma\eta}E\_{\mathrm{xc}} + \sum\_\sigma\int v\_{\mathrm{xc},\sigma}\left(\frac{\partial n\_\sigma}{\partial\epsilon\_{\gamma\eta}} + \frac12\frac{\partial n^{\mathrm{pcc}}}{\partial\epsilon\_{\gamma\eta}}\right) \\) | TN23 Eq. 46–48 |
| \\( E\_{\mathrm{xc}} \\)（GGA の追加分） | \\( -\Delta V\sum\_\sigma\sum\_p \frac{\partial f}{\partial\lVert\nabla n\_\sigma\rVert}\frac{\partial\_\gamma n\_\sigma\\, \partial\_\eta n\_\sigma}{\lVert\nabla n\_\sigma\rVert}\Big\rvert\_{\mathbf r\_p} \\) | TN23 付録 A |
| 格子モードの \\( E\_{\mathrm{na}} \\) | \\( \delta\_{\gamma\eta}E\_{\mathrm{na}} + \int V^{\mathrm{na}}\frac{\partial n}{\partial\epsilon\_{\gamma\eta}} + \int n \frac{\partial V^{\mathrm{na}}}{\partial\epsilon\_{\gamma\eta}} \\) | TN23 4.2 節と同じ考え方 |
| \\( U\_{\mathrm{core}} \\) | \\( -\frac12\sum\_{I,(J,\mathbf L)} Z\_I Z\_J\frac{\tau\_\gamma\tau\_\eta}{d^3} \\) | TN23 Eq. 57–59 |
| \\( U\_{\mathrm{H0}} \\) | \\( -\frac12\sum s\_{IJ}\frac{\partial J\_{IJ}}{\partial\tau\_\gamma}\tau\_\eta \\)（対称化） | TN23 Eq. 60–61 |
| `UpV` | \\( P\Omega\\,\delta\_{\gamma\eta} \\)（\\( P \\) は `MD.applied.pressure`） | 〔規約〕 |

- 格子上の関数のひずみ微分は、関数を元の格子へ引き戻して \\( \partial g((1+\boldsymbol\epsilon)(\mathbf r - \mathbf t))/\partial\epsilon\_{\gamma\eta} = \partial\_\gamma g\\,(\mathbf r - \mathbf t)\_\eta \\)
  とします（TN23 Eq. 36、38–40）。\\( \partial n/\partial\epsilon \\) には \\( \rho \\) の変化も入りますが、それは力と同じく Pulay の項にまとまります。
- \\( E\_{\delta\mathrm{ee}} \\) の最後の項は、\\( 1/\lvert\mathbf G\rvert^2 \\) のひずみ微分 \\( 2G\_\gamma G\_\eta/G^4 \\)（TN23 Eq. 43）から来ます。
  \\( \delta n(\mathbf G) = \Omega^{-1}\int\_\Omega \delta n\\, e^{-i\mathbf G\cdot\mathbf r} \\) です。
- GGA：既存実装は格子点の勾配を、格子ベクトル方向の中心差分 \\( \Delta\_i n = n(\mathbf r\_p + \mathbf h\_i) - n(\mathbf r\_p - \mathbf h\_i) \\)
  （\\( \mathbf h\_i = \mathbf a\_i/N\_i \\)）から、3×3 の連立方程式 \\( \mathbf h\_i\cdot\nabla n = \frac12\Delta\_i n \\) を解いて求めます（TN23 付録 A、Eq. 73–77）。
  この差分勾配では、格子点の値を固定してひずみをかけると勾配が \\( (1+\boldsymbol\epsilon)^{-T}\nabla n \\) に変わり、それが表の GGA の追加分を与えます。
  格子点の値の変化による部分は、差分勾配の転置を含む \\( v\_{\mathrm{xc}} \\) にまとまります（TN23 Eq. 91–93）。
- 交換相関の体積項は、原子密度の補正の原子中心の求積部分を含みません（原子中心の部分は対の相対位置だけに依存するので、
  他の対の項と同じく \\( \partial/\partial\tau\_\gamma\\,\tau\_\eta \\) で扱います）〔規約〕。

> **TN23 についての注記**：
> - Eq. 34 と Eq. 44 の体積項は \\( \delta\_{\gamma\eta}\int\delta n\\,\delta V\_H \\) と印字されていますが、
>   \\( E\_{\delta\mathrm{ee}} = \frac12\int\delta n\\,\delta V\_H \\) の積分領域の変化からは \\( \delta\_{\gamma\eta}E\_{\delta\mathrm{ee}} \\)（係数 \\( \frac12 \\) 付き）が出ます。
>   硬い Gauss 型電荷の周期配列で、ひずみをかけたエネルギーの有限差分と 3 項の和が \\( 10^{-10} \\) で一致するのは
>   \\( \frac12 \\) 付きの場合です（`tools/check_energy_force.py`）。既存実装も `UH1` を対角に 1 回加えます〔規約〕。
> - Eq. 57–59 の分子 \\( Z\_I - Z\_J \\) は \\( Z\_I Z\_J \\) の誤植と考えられます。また \\( \partial(1/\lvert\mathbf t\rvert)/\partial\epsilon\_{\gamma\eta} = -t\_\gamma t\_\eta/\lvert\mathbf t\rvert^3 \\) なので、
>   Eq. 59 の右辺には負号が必要です。上の表は負号付きです（既存実装も負号付き〔規約〕）。

### 出力と格子ベクトルの勾配〔規約〕

- `scf.stress.tensor on`（既定 off）、または格子を動かす最適化・NPT の MD で計算されます。
- 標準出力に「Stress tensor (Hartree/bohr^3)」として \\( \Sigma\_{\gamma\eta}/\Omega \\)（対称化後）の 3×3 行列が出ます。
  \\( \Sigma \\) は、単位胞が縮みたがる状態（引っ張られた状態）で対角成分が正になる符号です。
  圧力は \\( P = -\frac13\mathrm{tr}(\Sigma)/\Omega \\)（圧縮された状態で正）として使います。
- 軌道最適化、非共線スピン、DFT+U の 'full' 射影子とは併用できません（入力の検査で停止します）。
- 〔設計〕格子ベクトルを行とする行列 \\( A \\)（\\( A\_{ij} \\) は \\( \mathbf a\_i \\) の \\( j \\) 成分）では、ひずみで \\( A \to A(1+\boldsymbol\epsilon)^T \\) なので、
  格子ベクトルに関する勾配は

$$
\frac{\partial U\_{\mathrm{tot}}}{\partial A\_{ij}} = \sum\_k (A^{-1})\_{ki}\\, \Sigma\_{kj}
$$

  です。マニュアル 19.2 節が「格子ベクトルに関する勾配」と呼ぶ量に当たります。

〔未確定〕

- 出力ファイル `System.Name.out` に応力・格子ベクトルの勾配がどの書式で書かれるか。
- 出力される応力に `UpV` の寄与（加えた圧力）が含まれるかどうか（既定では圧力 0 なので影響しない）。
- 射影子展開を使わない組み合わせの \\( V^{\mathrm{na}} \\)、DFT+U、vdW の応力の詳細。

## テスト〔設計〕

### 入力の条件

- SCF の収束判定を十分に厳しくします（`scf.criterion` を \\( 10^{-10} \\) 程度）。力・応力は変分的でない誤差を持つため。
- 有限差分のテストは、絶縁体（例：H₂O 分子、Si 結晶）と低い電子温度で行うか、自由エネルギー \\( A \\) と比べます。
- 1 つの項だけを有効にして（他の項の寄与を 0 にして）比べるテストも用意します。既存実装の力の検査
  （`openmx -forcetest`、フラグ 0〜8 で運動・非局所・中性原子・差の Hartree・交換相関・電場・DFT+U を選ぶ、マニュアル 72 節）と同じ考え方です。

### 不変量と有限差分

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 有限差分 | 原子 1 個を各方向に \\( \pm 3\times 10^{-4} \\) bohr 動かした `Utot` の中心差分 \\( -[U(+)-U(-)]/(2h) \\) と解析的な力 | \\( 10^{-4} \\) Hartree/bohr（全項）。項ごとでは \\( 10^{-5} \\) |
| 有限差分 | 二中心・三中心の項だけ（`Ukin`、`Unl`、射影子展開の `Una`）と Pulay の項：固定した \\( \rho \\)、\\( \varepsilon \\) での微分 | \\( 10^{-7} \\) |
| 有限差分 | `Ucore + UH0` だけ | \\( 10^{-6} \\)（円筒格子の求積誤差の程度） |
| 有限差分 | ひずみ \\( \pm 10^{-4} \\) の `Utot` の中心差分と \\( \Sigma\_{\gamma\eta} \\)（6 成分）。非対角は \\( \epsilon\_{\gamma\eta} = \epsilon\_{\eta\gamma} = h/2 \\) で | \\( 10^{-4} \\) Hartree（単位胞あたり） |
| 有限差分 | 一様な体積変化 \\( E(V) \\) の微分と \\( -P \\) | 上と同程度 |
| 不変量 | 力の総和：二中心・三中心・Pulay・`Ucore + UH0` の寄与はそれぞれ総和 0 | \\( 10^{-10} \\) |
| 不変量 | 力の総和：格子上の寄与の総和（格子が原子と一緒に動かないための誤差） | 格子のカットオフとともに減少することを確認。150 Ry で \\( 10^{-3} \\) Hartree/bohr 程度を目安 |
| 不変量 | 全原子を格子点の間隔の整数倍だけ平行移動：エネルギー・力が不変 | \\( 10^{-10} \\) |
| 不変量 | 全原子を任意の量だけ平行移動：エネルギーの変化（egg-box 効果）が格子のカットオフとともに減少 | 傾向を確認 |
| 不変量 | 孤立分子の回転：エネルギーは egg-box 効果の範囲で不変、力は同じ回転で変換 | egg-box 効果の範囲 |
| 不変量 | 応力テンソルが対称 | 丸め誤差 |
| 恒等式 | \\( U\_{\mathrm{ele}} = g\_s\sum\rho H \\) | \\( 10^{-10} \\) |
| 恒等式 | 収束後の Harris 型の `Utot` と通常の `Utot` | \\( 10^{-8} \\) |
| 恒等式 | 孤立原子：`Ucore` = 0、`UH0` = 自己項 | \\( 10^{-10} \\) |
| 恒等式 | 原子対の範囲の端で、遮蔽された対のエネルギーが 0 | \\( 10^{-10} \\) |
| 解析 | Gauss 型の密度での \\( J\_{IJ}(d) \\) と \\( Z\_I Z\_J\\,\mathrm{erf}(\sqrt\mu d)/d \\) | \\( 10^{-8} \\)（相対） |
| 解析 | 一般化固有値問題で \\( \partial\sum fe = \mathrm{Tr}[\rho\partial H] - \mathrm{Tr}[\varepsilon\partial S] \\) | \\( 10^{-8} \\) |
| 解析 | 周期的な硬い Gauss 型電荷の Hartree 応力 | \\( 10^{-8} \\) |

最後の 3 つの解析テストは `tools/check_energy_force.py` に実装しています。

### fixture との照合

既存実装から、次の量を書き出して比較します。

- エネルギー：`Uele`、`Ukin`、`UH0`、`UH1`、`Una`、`Unl`、`Uxc0`、`Uxc1`、`Ucore`、`UvdW`、`Utot`（出力ファイルの値）。
  `scf.Exc0.correction.flag` の on・off、`scf.ProExpn.VNA` の on・off の両方。
- 力：`<coordinates.forces` の値と、項ごとの寄与（`-forcetest` の各フラグでの解析的な力）。
- 応力：標準出力の応力テンソル。
- 中間量：収束後の \\( \rho \\)、\\( \varepsilon \\)（原子対ごとのブロック）、格子上の \\( \delta V\_H \\)、\\( v\_{\mathrm{xc},\sigma} \\)、\\( n^{a} \\)。
  これらを使うと、固有値問題を解かずにエネルギーの各成分と力を再計算して比べられます。
- 対象：H₂O 分子（`Cluster`）、Si 結晶（`Band`、応力を含む）、bcc Fe（共線スピン、`Band`）。

| 量 | 許容誤差の目安（同じ格子・同じ基底のとき） |
|---|---|
| 各エネルギー成分 | \\( 10^{-6} \\) Hartree |
| `Utot` | \\( 10^{-6} \\) Hartree |
| 力 | \\( 10^{-5} \\) Hartree/bohr |
| 応力 \\( \Sigma/\Omega \\) | \\( 10^{-7} \\) Hartree/bohr³ |

格子点の数や原点の取り方が既存実装と異なる場合は、エネルギーの差は egg-box 効果の程度まで広がるので、
許容誤差は [グリッド](grid-poisson-xc.md) の章の格子の規則が確定してから決め直します。

## 〔未確定〕

- 非共線スピン・スピン軌道相互作用でのエネルギー成分と力（\\( \rho \\) の複素数部分の寄与）。
- `Uhub`、`Ucs`、`Uzs`、`Uzo`、`Uch` の定義と力・応力。
- 「開いた内殻」の扱いを持つ擬ポテンシャルでの部分内殻補正の力（スピンの片方の \\( v\_{\mathrm{xc}} \\) だけを使う場合がある）。
- GGA の \\( v\_{\mathrm{xc}} \\) が差分勾配の転置を含む離散的に厳密な微分になっているか。
- 原子密度の補正で、原子中心の求積に部分内殻補正の密度を含めないことを再実装でも踏襲するか。
- 中性原子項のうち射影子展開を使わない組み合わせの行列要素の微分。
- 応力・格子ベクトルの勾配の、出力ファイル上の書式。
