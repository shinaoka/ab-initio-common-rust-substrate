# 解析機能：バンド・状態密度・電荷解析・交換相互作用・Berry 位相・光学応答

> 状態：初版（レビュー待ち）。式と手法は、マニュアル 4.0（22・23・25・33–35・44–45・50・53・56・58・59・62・70–71 節）、
> 技術ノート（TN9「Charge Analysis」、TN11「Electric Polarization by Berry Phase」、「Calculation of Optical Matrix Elements」、
> 「Dielectric Function and Conductivity」、「Analysis of spin texture in the k-space」）と論文
> （Han–Ozaki–Yu, PRB **70**, 184421 (2004)；Terasawa ら, JPSJ **88**, 114706 (2019)；Sawahata ら, JJAP **57**, 030309 (2018)；
> Kotaka ら, JJAP **52**, 035204 (2013)）を、このノートの記法で書き直したものです。
> ファイルの配置、既定値、数値上の選択は既存実装の挙動として確認したもので、〔規約〕と記します。
> Mulliken の和則、四面体法、Berry 位相、Chern 数、交換結合の式の関係は `tools/check_analysis.py` で数値的に確認しています。

## この章の範囲

SCF が終わった後の解析（後処理）を扱います。優先度の高い順に次のとおりです。

1. 後処理の入力となる `.scfout` ファイルの形式と、OpenMX API（1 節）。
2. 電荷解析：Mulliken・Voronoi・ESP 電荷（4 節）。
3. 状態密度（DOS）と部分状態密度（PDOS）、後処理プログラム DosMain（3 節）。
4. 交換結合定数 \\( J\_{ij} \\)（jx、6 節）。
5. Berry 位相による電気分極（polB、7 節）。

このほか、バンド構造（2 節）、軌道分解全エネルギー（5 節）、Chern 数・Berry 曲率・\\( \mathbb Z\_2 \\) 不変量（8 節）、
k 空間のスピンテクスチャ（9 節）、光学伝導度と誘電関数（10 節）、格子上の量の出力（cube ファイル、分子軌道、STM、仕事関数、11 節）を簡潔に扱います。

非共線計算の Mulliken のスピンモーメントと軌道モーメントの定義は [非共線 DFT のエネルギーと力](noncollinear-energy-force.md) の 6 節にあり、ここでは繰り返しません。
NEGF の出力（透過率、固有チャネルなど）は [NEGF](negf.md) にあります。自然軌道解析（NPA、NBO、LNAO）、CWF・MLWF、バンドのアンフォールディング、
COHP は対象外とします（〔未確定〕の節を参照）。

## 記法〔設計〕

[k 点と Bloch 和](kpoints-bloch.md) と [密度行列と電荷ミキシング](density-mixing.md) の記法に従います。

| 記号 | 意味 |
|---|---|
| \\( O\_{I\alpha,J\beta}(\mathbf L) \\) | 実空間行列（[疎ブロック形式](../design/sparse-block-format.md)）。\\( \mathbf L \\) は列側の原子 \\( J \\) の周期像 |
| \\( O(\mathbf k) = \sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L} O(\mathbf L) \\) | Bloch 和（格子ベクトルだけの位相） |
| \\( c\_{a\mu}(\mathbf k) \\) | 一般化固有値問題の固有ベクトル（\\( c^\dagger S(\mathbf k) c = 1 \\)）。\\( a = (I\alpha) \\) は全軌道の番号 |
| \\( \rho^{\sigma}(\mathbf L) \\) | 実空間の密度行列（非スピン偏極では 1 スピン分） |
| \\( \mu \\)、\\( T \\) | 化学ポテンシャル（Hartree）、電子温度（K） |
| \\( f(\varepsilon) \\) | Fermi 分布 \\( 1/(1+e^{(\varepsilon-\mu)/k\_{\mathrm B}T}) \\) |
| \\( g\_s \\) | スピンの縮重度（非スピン偏極で 2、それ以外で 1） |
| \\( \Omega \\) | 単位胞の体積（bohr³） |
| \\( \mathbf a\_i \\)、\\( \mathbf b\_i \\) | 格子ベクトルと逆格子ベクトル（\\( \mathbf a\_i\cdot\mathbf b\_j = 2\pi\delta\_{ij} \\)） |
| \\( D^{d}\_{I\alpha,J\beta}(\mathbf L) \\) | 位置演算子の行列 \\( \langle\phi\_{I\alpha}\rvert (r\_d - R\_{I,d}) \lvert\phi\_{J\beta}(\mathbf L)\rangle \\)（\\( d = x, y, z \\)、原子 \\( I \\) の中心から測る） |
| \\( \nabla^{d}\_{I\alpha,J\beta}(\mathbf L) \\) | 微分演算子の行列 \\( \langle\phi\_{I\alpha}\rvert \partial\_d \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) |

\\( \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) は \\( \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L) \\) の略です。

## 1. `.scfout` ファイルと OpenMX API

### 1.1 役割〔公開〕

マニュアル 70.1 節：入力キーワード `HS.fileout`（`on`／`off`、既定 `off`）を `on` にすると、計算の最後に
`System.Name.scfout` が書き出されます。Hamiltonian、重なり行列、密度行列、位置演算子と運動量演算子の行列要素を含み、
後処理プログラム（jx、polB、kSpin、Z2FH、calB、および読み取りの例）はこのファイルだけを入力とします。
マニュアルはファイルの配置を示さず、読み取りの例のソースを参照するよう述べるだけです。以下の配置は既存実装から確認した事実です。

### 1.2 バイナリーの配置〔規約〕

- 書式なしのバイナリーで、整数は 4 バイト、実数は 8 バイト（倍精度）、書き手の計算機のバイト順です。レコード区切りはありません。
- 以下の「ブロック行列」は、原子 \\( I = 1, \dots, N\_{\mathrm{atom}} \\) の順に、その近接リストの各要素 \\( h = 0, \dots, n\_I \\) について、
  \\( n\_{\mathrm{orb}}(I)\times n\_{\mathrm{orb}}(J\_h) \\) のブロックを **行優先**（行 = 原子 \\( I \\) の軌道、列 = 相手 \\( J\_h \\) の軌道。列が速く変わる）で並べたものです。
  1 つのブロック行列の総要素数は \\( \sum\_I \sum\_{h=0}^{n\_I} n\_{\mathrm{orb}}(I)\\, n\_{\mathrm{orb}}(J\_h) \\) です。

| 順 | 内容 | 型 × 個数 |
|---|---|---|
| 1 | 見出し：\\( N\_{\mathrm{atom}} \\)、\\( s + 4v \\)、NEGF 用の原子数 3 つ、\\( N\_{\mathrm{cell}} \\) | int × 6 |
| 2 | 位置演算子の次数の最大値 \\( p\_{\max} \\)（OpenMX 4.0 では 1） | int × 1 |
| 3 | 周期像の Cartesian 並進ベクトル（bohr）、\\( r = 0, \dots, N\_{\mathrm{cell}} \\) の各 \\( r \\) について \\( (\cdot, x, y, z) \\) | double × \\( 4(N\_{\mathrm{cell}}+1) \\) |
| 4 | 周期像の整数座標、各 \\( r \\) について \\( (\cdot, n\_1, n\_2, n\_3) \\) | int × \\( 4(N\_{\mathrm{cell}}+1) \\) |
| 5 | 各原子の軌道数 \\( n\_{\mathrm{orb}}(I) \\) | int × \\( N\_{\mathrm{atom}} \\) |
| 6 | 各原子の近接原子数 \\( n\_I \\)（自分自身を除く） | int × \\( N\_{\mathrm{atom}} \\) |
| 7 | 近接原子の番号（1 始まり）：原子ごとに \\( n\_I + 1 \\) 個 | int × \\( \sum\_I (n\_I+1) \\) |
| 8 | 近接原子の周期像の番号 \\( r \\)：原子ごとに \\( n\_I + 1 \\) 個 | int × \\( \sum\_I (n\_I+1) \\) |
| 9 | 格子ベクトル \\( \mathbf a\_1, \mathbf a\_2, \mathbf a\_3 \\)（bohr）、それぞれ \\( (\cdot, x, y, z) \\) | double × 12 |
| 10 | 逆格子ベクトル \\( \mathbf b\_1, \mathbf b\_2, \mathbf b\_3 \\)（\\( 2\pi \\) を含む、bohr⁻¹）、それぞれ \\( (\cdot, x, y, z) \\) | double × 12 |
| 11 | 原子位置（Cartesian、bohr）、原子ごとに \\( (\cdot, x, y, z) \\) | double × \\( 4N\_{\mathrm{atom}} \\) |
| 12 | Hamiltonian \\( H \\)：スピン成分 \\( \sigma = 0, \dots, s \\) ごとにブロック行列 1 つ（Hartree） | double |
| 13 | （\\( s = 3 \\) のときだけ）Hamiltonian の虚部の追加分：3 成分のブロック行列 | double |
| 14 | 重なり行列 \\( S \\)：ブロック行列 1 つ | double |
| 15 | 位置演算子 \\( D^{d} \\)：\\( d = x, y, z \\) の順、各 \\( d \\) で次数 \\( 1, \dots, p\_{\max} \\) の順にブロック行列 | double |
| 16 | 微分演算子 \\( \nabla^{d} \\)：\\( d = x, y, z \\) の順にブロック行列 | double |
| 17 | 密度行列 \\( \rho \\)：スピン成分 \\( \sigma = 0, \dots, s \\) ごとにブロック行列 1 つ | double |
| 18 | 密度行列の虚部の追加分：**常に** 2 成分のブロック行列（共線計算でも書かれる） | double |
| 19 | 解法の番号（下表） | int × 1 |
| 20 | \\( \mu \\)（Hartree）、\\( T \\)（K）、芯電荷の双極子 3 成分（Debye）、一様背景電荷の双極子 3 成分（Debye）、価電子数 \\( N\_{\mathrm e} \\)、全スピン \\( S \\) | double × 10 |
| 21 | 入力ファイルの行数 \\( n\_{\mathrm{line}} \\) | int × 1 |
| 22 | 入力ファイルの各行：256 バイトの固定長レコード \\( n\_{\mathrm{line}} \\) 個 | char × \\( 256\\,n\_{\mathrm{line}} \\) |

各項目の意味：

- **見出しの第 2 要素**：スピンの種類 \\( s \\)（0：非スピン偏極、1：共線スピン偏極、3：非共線）と形式の版 \\( v \\) を \\( s + 4v \\) にまとめた値です。
  OpenMX 3.9・4.0 は \\( v = 3 \\) で、値は 12、13、15 のいずれかです。読み手は \\( s = (\cdot) \bmod 4 \\)、\\( v = \lfloor(\cdot)/4\rfloor \\) と分解し、
  \\( v \ne 3 \\) なら読み取りを拒否します（\\( v = 0 \\) は 3.8 以前、1・2 は開発版の形式）。
- **バイト順の判定**：既存の読み手は、第 2 要素が \\( 4\times 3 + 3 = 15 \\) を超えるとき、ファイル全体をバイト順の逆転したものとみなして各値を反転して読みます。
  反転後も範囲外なら誤りとします。
- **NEGF 用の原子数**：中央・左・右の領域の原子数で、NEGF 以外の計算では 0 です。読み手は使いません。
- **周期像**：\\( N\_{\mathrm{cell}} = (2N\_c+1)^3 - 1 \\) で、\\( N\_c \\) は近接原子をすべて含むのに十分な各方向の最大の整数座標です。
  番号 \\( r = 0 \\) が原点の胞、\\( r \ge 1 \\) は \\( n\_1 \\) を最も外側、\\( n\_3 \\) を最も内側として \\( [-N\_c, N\_c]^3 \\) から原点を除いて辞書式に並べたものです。
  読み手は番号の並びに依存せず、項目 4 の整数座標を使うべきです。
- **4 要素のベクトル**：項目 3・4・9・10・11 は先頭の要素を使わず（内容は不定）、要素 1–3 が \\( x, y, z \\)（または \\( n\_1, n\_2, n\_3 \\)）です。
- **近接リスト**：原子 \\( I \\) の第 0 要素は自分自身（原子番号 \\( I \\)、周期像 \\( r = 0 \\)）です。第 1 要素以降は
  \\( \lVert \mathbf R\_J + \mathbf L - \mathbf R\_I\rVert \le r\_{c,I} + r\_{c,J} \\)（等号を含む）を満たすすべての \\( (J, \mathbf L) \\) で、並び順は規定しません。
  \\( (J, \mathbf L) \\) が \\( I \\) の近接なら \\( (I, -\mathbf L) \\) は \\( J \\) の近接です。これは [疎ブロック形式](../design/sparse-block-format.md) の原子対の集合と同じです。
- **項目 12 の \\( H \\)**：共線計算では \\( \sigma = 0 \\)（上向き）、\\( \sigma = 1 \\)（下向き、\\( s = 1 \\) のみ）の実数行列です。
  非共線計算（\\( s = 3 \\)）では 4 成分 \\( (H\_0, H\_1, H\_2, H\_3) \\) で、項目 13 の 3 成分 \\( (X^{\uparrow\uparrow}, X^{\downarrow\downarrow}, X^{\uparrow\downarrow}) \\)
  （スピン軌道相互作用、DFT+U、スピン拘束、Zeeman 項などによる虚部）と合わせて

  $$
  H^{\uparrow\uparrow} = H\_0 + iX^{\uparrow\uparrow},\qquad H^{\downarrow\downarrow} = H\_1 + iX^{\downarrow\downarrow},\qquad
  H^{\uparrow\downarrow} = H\_2 + i\left(H\_3 + X^{\uparrow\downarrow}\right),\qquad H^{\downarrow\uparrow}\_{J\beta,I\alpha}(-\mathbf L) = \overline{H^{\uparrow\downarrow}\_{I\alpha,J\beta}(\mathbf L)}
  $$

  となります。k 空間の \\( 2N\times 2N \\) 行列の組み立ては [k 点と Bloch 和](kpoints-bloch.md) のとおりです。
- **項目 12 と 17 の時点**：計算（構造最適化・MD を含む）の最後に書かれます。\\( H \\) は最後に対角化した Hamiltonian（最後の SCF 反復の入力密度から作ったもの）、
  \\( \rho \\) はそれを対角化して得た出力の密度行列です。SCF が収束していれば両者は収束判定の精度で整合します。
- **軌道の縮約**：軌道最適化（縮約）を使った場合、\\( n\_{\mathrm{orb}} \\) と \\( H \\)、\\( S \\) は縮約後の基底のものです。
- **項目 15 の \\( D \\)**：\\( \langle\phi\_{I\alpha}\rvert (r\_d - R\_{I,d})^{p} \lvert\phi\_{J\beta}(\mathbf L)\rangle \\)（\\( p = 1, \dots, p\_{\max} \\)）で、原点は **行側の原子** \\( I \\) の中心です。
  解析的ではなく、Poisson 方程式と同じ実空間格子上の和（格子点の体積を掛けたもの）で評価した値で、精度は格子の細かさ（`scf.energycutoff`）に依存します。
- **項目 16 の \\( \nabla \\)**：\\( \langle\phi\_{I\alpha}\rvert \partial\_d \lvert\phi\_{J\beta}(\mathbf L)\rangle \\) です（運動量演算子の行列はこれの \\( -i \\) 倍）。これも格子上の和です。
- **項目 17・18 の \\( \rho \\)**：共線計算では [密度行列と電荷ミキシング](density-mixing.md) の \\( \rho^\sigma(\mathbf L) \\)（非スピン偏極では 1 スピン分）で、項目 18 は 0 です。
  非共線計算では項目 17 が \\( (\mathrm{Re}\\,\rho^{\uparrow\uparrow}, \mathrm{Re}\\,\rho^{\downarrow\downarrow}, \mathrm{Re}\\,\rho^{\uparrow\downarrow}, \mathrm{Im}\\,\rho^{\uparrow\downarrow}) \\)、
  項目 18 が \\( (\mathrm{Im}\\,\rho^{\uparrow\uparrow}, \mathrm{Im}\\,\rho^{\downarrow\downarrow}) \\) です。
- **項目 19 の解法の番号**：`scf.EigenvalueSolver` の値に対応します。

  | 値 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | 解法 | Recursion | Cluster | Band | NEGF | DC | GDC | Cluster-DIIS | Krylov | Cluster2 | EGAC | DC-LNO | Cluster-LNO |

  `Band` で `scf.Kgrid` が `1 1 1` のとき、計算は Cluster の経路（周期境界つき Γ 点）で行われますが（[k 点と Bloch 和](kpoints-bloch.md)）、ファイルには 3（Band）が書かれます。
- **項目 20**：\\( T \\) は `scf.ElectronicTemperature` の値（K）です。双極子は Debye 単位（1 \\( e \\)·bohr = 2.54174776 Debye）で、
  芯電荷の項は \\( \sum\_I Z\_I^{\mathrm{val}} \mathbf R\_I \\)、背景電荷の項は電荷ドープ（`scf.system.charge` \\( = Q \\)）の一様背景 \\( Q/\Omega \\) についての
  \\( -\int\_\Omega (Q/\Omega)\\,\mathbf r\\,d^3r \\) を格子上の和で評価したものです。全スピン \\( S \\) は共線で Mulliken の \\( \frac12\sum\_I (M^\uparrow\_I - M^\downarrow\_I) \\) です（4.1 節）。
- **項目 22**：各行は、入力ファイルを 255 文字までの行単位で読んだ文字列（改行文字を含み、終端の NUL 以降の内容は不定）です。
  255 文字を超える行は複数のレコードに分かれます。読み手は NUL までを取り出して連結し、入力ファイルを復元します。
  入力ファイル名が与えられない場合は、項目 21・22 がありません。

### 1.3 疎ブロック形式との対応〔設計〕

`.scfout` の内容は、このノートの [疎ブロック形式](../design/sparse-block-format.md) にそのまま写せます。

- 原子対 \\( (I, J\_h, \mathbf n\_h) \\)：\\( J\_h \\) は項目 7、\\( \mathbf n\_h \\) は項目 8 の \\( r \\) を項目 4 で整数座標に直したもの。
- ブロックの値：項目 12–18 の各ブロック（行優先）。疎ブロック形式のメモリー配置が列優先なら、読み込み時に転置します。
- 非共線：上の式で複素数のスピンブロックにまとめます。

新しい実装は `.scfout` の **読み手と書き手** を、疎ブロック形式との変換器として持ちます（以下の解析はすべて疎ブロック形式を入力とする）。
書き手は既存の後処理プログラムとの互換のために、版 3 の配置を厳密に再現します。

### 1.4 OpenMX API〔公開〕〔規約〕

〔公開〕マニュアル 71 節：OpenMX を外部のプログラム（C、Fortran、Python）からライブラリーとして呼ぶ仕組みで、
静的ライブラリー（C・Fortran）と Python の拡張モジュールが用意されています。典型的な使い方は、MPI のコミュニケーターを分割して、
複数の入力ファイルを並列に計算することです。

〔規約〕公開されている入口は 1 つで、引数は (MPI コミュニケーター、入力ファイル名、OpenMP のスレッド数) の 3 つ、戻り値は整数（0 が正常終了）です。
Fortran 版はコミュニケーターを Fortran のハンドルで受け取り、Python 版は同じ 3 つの引数をとる関数 1 つです。

- 呼ぶと、渡したコミュニケーターの上で、入力ファイルに書かれた計算（SCF、構造最適化、MD、NEB など）を通常の実行と同じく行い、
  出力はすべて通常のファイル（`.out`、`.scfout`、cube ファイルなど）に書かれます。
- メモリー上の行列や物理量を呼び出し側に返す関数はありません。後処理に使う量は `.scfout` などのファイルを介して受け取ります。

〔設計〕新しい実装の API は、ファイルを介さずに疎ブロック形式の \\( H \\)、\\( S \\)、\\( \rho \\) と構造を返す関数を持ちます。
既存の API と同じ「入力ファイルを渡して実行する」入口は、互換のための薄い層として用意します。

## 2. バンド構造

### 2.1 入力〔公開〕〔規約〕

〔公開〕マニュアル 22 節：`Band.dispersion`（`on`／`off`、既定 `off`）、経路の数 `Band.Nkpath`、
経路のブロック `<Band.kpath`（各行：点の数 \\( N \\)、始点 \\( (q\_1, q\_2, q\_3) \\)、終点 \\( (q\_1, q\_2, q\_3) \\)、始点と終点の記号）で指定します。
`<Band.KPath.UnitCell` で、k 点の座標を解釈するための別の単位胞（面心・体心格子の慣用胞など）を与えられます。省略時は `Atoms.UnitVectors` です。

〔規約〕

- 座標は逆格子ベクトルに対する分数座標です。`<Band.KPath.UnitCell` を与えた場合は、その胞の逆格子ベクトルで書かれた座標を、
  計算に使う胞の分数座標に変換します（\\( \mathbf k \\) を Cartesian に直してから \\( q\_i = \mathbf k\cdot\mathbf a\_i/2\pi \\)）。
- 各経路の点は \\( \mathbf q\_j = \mathbf q\_{\mathrm s} + (\mathbf q\_{\mathrm e} - \mathbf q\_{\mathrm s})\\,(j-1)/(N-1) \\)、\\( j = 1, \dots, N \\) で、
  **両端を含み**、分数座標で等間隔です。隣り合う経路の共通の端点は 2 回現れます。\\( N = 1 \\) は不可です（0 除算）。
- 各点で [一般化固有値問題](generalized-eigen.md) を解き、全固有値を出します（後処理の規則 A′）。

### 2.2 出力 `System.Name.Band`〔規約〕

テキストで、次の順に書かれます。

1. \\( n\\;\\; s\\;\\; \mu \\)：\\( n \\) は全軌道数（非共線では \\( 2N \\) と書き、\\( s \\) は 0 と書く）、\\( \mu \\) は Hartree。
2. 逆格子ベクトル 9 個（\\( \mathbf b\_1, \mathbf b\_2, \mathbf b\_3 \\) の Cartesian 成分、\\( 2\pi \\) を含む、bohr⁻¹）。
3. 経路の数。続いて経路ごとに \\( N \\)、始点、終点、記号 2 つ。
4. 経路・点・スピンの順に、\\( n\\;\\; q\_1\\;\\; q\_2\\;\\; q\_3 \\) の行と、\\( n \\) 個の固有値（Hartree、\\( \mu \\) を引かない）の行。

〔公開〕変換プログラム bandgnu13 は、横軸を経路に沿った距離、縦軸を \\( \varepsilon - \mu \\)（eV）にした gnuplot 用のデータを作ります（マニュアル 22 節）。

## 3. 状態密度と部分状態密度

### 3.1 計算の段階〔公開〕〔規約〕

〔公開〕マニュアル 23.1 節：`Dos.fileout on` で、SCF の後に DOS 用の固有値と固有ベクトルの情報を `System.Name.Dos.val`（テキスト）と
`System.Name.Dos.vec`（バイナリー）に書き出し、別のプログラム DosMain で DOS・PDOS を求めます。
`Dos.Erange` はエネルギーの範囲（eV、\\( \mu \\) を原点とする）、`Dos.Kgrid` は k 点の分割数です。O(N) 法でも DOS を出せます（Gauss 型の広がりのみ）。
23.2 節：`DosGauss.fileout on` は、固有ベクトルを書かずに Gauss 型の広がりで DOS をその場で求める方式で、`Dos.fileout` と同時には使えません。

〔規約〕

| キーワード | 既定値 | 備考 |
|---|---|---|
| `Dos.fileout` | `off` | |
| `Dos.Erange` | `-20.0 20.0`（eV） | \\( \mu \\) からの相対値 |
| `Dos.Kgrid` | `scf.Kgrid` の値 | |
| `DosGauss.fileout` | `off` | Band（周期的な Γ 点計算を含む）でのみ使える |
| `DosGauss.Num.Mesh` | 200 | エネルギーの分割数 |
| `DosGauss.Width` | 0.2（eV） | Gauss 関数の幅 \\( a \\) |

Band の場合の DOS 用の計算：

- k 点は \\( q\_{i,l} = l/N\_i \\)（\\( l = 0, \dots, N\_i - 1 \\)）の Γ 中心の格子で、`scf.Generation.Kpoint` によらずこの形です
  （分割数 2 以上の方向に [k 点と Bloch 和](kpoints-bloch.md) と同じ微小なずれ \\( \delta \\) が加わります）。時間反転による削減は **しません**。
  四面体法で周期的に隣の点を参照するため、格子全体が必要です。
- 各 k 点・スピンで固有値問題を解き（正則化は規則 A′）、エネルギーの窓に入る状態の番号の範囲 \\( [l\_{\min}, l\_{\max}] \\) を全 k 点で共通に決めます。
  各 k 点で \\( \mu + E\_{\min} \\) を超える最初の状態の 1 つ前から、\\( \mu + E\_{\max} \\) を超える最初の状態までをとり、全 k 点で和集合をとった後、
  両側に \\( \max(\lfloor n/20\rfloor, 10) \\) 個の余裕を加えます（\\( 1 \le l \le n \\) に切り詰める）。
- 状態 \\( \mu \\)（k 点 \\( \mathbf k \\)、スピン \\( \sigma \\)）の各軌道への重みは Mulliken の総和型の重み

  $$
  w\_{a\mu}(\mathbf k) = \mathrm{Re}\left[\overline{c\_{a\mu}(\mathbf k)}\\,\bigl(S(\mathbf k)\\,c\_{\mu}(\mathbf k)\bigr)\_a\right],\qquad \sum\_a w\_{a\mu}(\mathbf k) = 1
  $$

  です（4.1 節の Mulliken 電荷と同じ分け方）。

### 3.2 `.Dos.val` と `.Dos.vec` の形式〔規約〕

`.Dos.val` はキーワード形式のテキストで、次の項目を持ちます。

| キーワード | 内容 |
|---|---|
| `mode` | 1：Cluster・Band（固有状態ごと）、5：DC、6：NEGF、7：`DosGauss.fileout`（エネルギー格子ごと） |
| `NonCol` | 0：共線、1：非共線 |
| `N` | 全軌道数 \\( n \\) |
| `Nspin` | \\( s \\)（非共線では 1 と書く） |
| `Erange` | `Dos.Erange`（Hartree） |
| `irange` | \\( l\_{\min}\\; l\_{\max} \\)（1 始まり） |
| `Kgrid` | \\( N\_1\\; N\_2\\; N\_3 \\)（Cluster では `1 1 1`） |
| `atomnum` | 原子数 |
| `<WhatSpecies` … `WhatSpecies>` | 各原子の種の番号（0 始まり） |
| `SpeciesNum` | 種の数 |
| `<Spe_Total_CNO` … `Spe_Total_CNO>` | 各種の軌道数 |
| `MaxL` | 角運動量の最大値 |
| `<Spe_Num_CBasis` … `Spe_Num_CBasis>` | 各種・各 \\( l \\) の動径関数（多重度）の数 |
| `ChemP` | \\( \mu \\)（Hartree） |
| `<SpinAngle` … `SpinAngle>` | （非共線のみ）各原子の Mulliken スピンの向き \\( (\theta, \varphi) \\)（ラジアン） |
| `<Eigenvalues` … `Eigenvalues>` | k 点・スピンごとに 1 行：\\( l\_1\\; l\_2\\; l\_3 \\)（k 点の整数番号）に続いて \\( \varepsilon\_{l\_{\min}}, \dots, \varepsilon\_{l\_{\max}} \\)（Hartree） |

`.Dos.vec` は、k 点・スピン・状態 \\( l\_{\min} \dots l\_{\max} \\) の順のレコードの列で、各レコードは int × 3（k 点の整数番号）と
4 バイト実数 × \\( n \\)（重み \\( w\_{a\mu} \\)、軌道の順）です。k 点の並びは計算の分担で決まるので、読み手は整数番号で k 点を識別します。
非共線では 4 バイト実数 × \\( 2n \\) で、原子ごとに、その原子の Mulliken スピンの向きを量子化軸とした「多数スピン側」の重み \\( n\_{\mathrm{orb}}(I) \\) 個と
「少数スピン側」の重み \\( n\_{\mathrm{orb}}(I) \\) 個が並びます。

### 3.3 DosMain の計算〔公開〕〔規約〕

〔公開〕マニュアル 23.1 節：DosMain は `.Dos.val` と `.Dos.vec` を読み、方法（四面体法または Gauss 型の広がり）、
全 DOS か PDOS か、PDOS の対象原子、Gauss 型の幅 \\( a \\)（eV、関数 \\( \exp(-(E/a)^2) \\)）を対話的に尋ねます。
出力は 1 列目がエネルギー（eV、\\( \mu \\) が原点）、2 列目が DOS（eV⁻¹）、3 列目が積分値です。スピン偏極では上向き・下向きの DOS と積分値が並びます。

〔規約〕

- **エネルギー格子**：\\( E\_j = E\_{\min} + (E\_{\max} - E\_{\min})\\,j/(N\_E - 1) \\)、\\( j = 0, \dots, N\_E - 1 \\)（\\( \mu \\) 基準）。
  四面体法では \\( N\_E = \lfloor (E\_{\max} - E\_{\min})/1\\,\mathrm{meV}\rfloor \\)、Gauss 型では \\( N\_E = \lfloor 5(E\_{\max} - E\_{\min})/a\rfloor \\) です。
- **Gauss 型**：

  $$
  D\_\sigma(E) = \frac{1}{N\_k}\sum\_{\mathbf k}\sum\_{\nu} W\_{\nu\sigma}(\mathbf k)\\, \frac{1}{a\sqrt\pi}\exp\left[-\left(\frac{E + \mu - \varepsilon\_{\nu\sigma}(\mathbf k)}{a}\right)^2\right]
  $$

  で、全 DOS では \\( W = \sum\_a w\_a \\)（= 1）、PDOS では対象の軌道の重みです。Gauss 関数は中心から \\( \pm 3a \\) 程度（格子点で数えて \\( 3a/\Delta E + 3 \\) 点）で打ち切ります。
  単位は 1 スピン・1 単位胞あたり eV⁻¹ です（\\( \nu \\) は状態の番号、\\( \mu \\) は化学ポテンシャル）。
- **四面体法**：k 点の格子の各小立方体（8 頂点、頂点の番号を \\( 4i + 2j + k \\)、\\( i, j, k \in \\{0, 1\\} \\) とする）を、対角 0–7 を共有する 6 個の四面体
  \\( \\{0,1,2,5\\}, \\{1,2,3,5\\}, \\{2,3,5,7\\}, \\{0,2,4,5\\}, \\{2,4,5,6\\}, \\{2,5,6,7\\} \\) に分けます（周期的に折り返す）。
  各四面体の頂点のエネルギーを昇順 \\( e\_1 \le e\_2 \le e\_3 \le e\_4 \\) にして、線形補間の DOS（Blöchl ら, PRB **49**, 16223 (1994) の Appendix C、四面体あたり 1 に規格化）

  $$
  g(E) = \begin{cases}
  \dfrac{3(E - e\_1)^2}{e\_{21}e\_{31}e\_{41}} & e\_1 < E < e\_2 \\\\[2mm]
  \dfrac{1}{e\_{31}e\_{41}}\left[3e\_{21} + 6(E - e\_2) - \dfrac{3(e\_{31} + e\_{42})(E - e\_2)^2}{e\_{32}e\_{42}}\right] & e\_2 < E < e\_3 \\\\[2mm]
  \dfrac{3(e\_4 - E)^2}{e\_{41}e\_{42}e\_{43}} & e\_3 < E < e\_4
  \end{cases}
  $$

  （\\( e\_{ij} = e\_i - e\_j \\)、範囲外は 0）を足し、\\( 6N\_k \\) で割ります。Blöchl の補正項は使いません。
  全 DOS では各状態の重みを小立方体の基準の頂点の値で一定とし、PDOS では頂点の重みを四面体内で線形補間した式（同 Appendix B を微分したもの）を使います。
- **積分値**：DOS を Simpson 則で \\( E\_0 \\) から \\( E\_j \\) まで積分した値です。
- **出力**：非スピン偏極では \\( 2D \\) と積分値を出します。スピン偏極では \\( D\_\uparrow \\)、\\( -D\_\downarrow \\)（**負の符号を付ける**）、\\( \int D\_\uparrow \\)、\\( \int D\_\downarrow \\) の順です。
  ファイル名は `System.Name.DOS.Tetrahedron`、`System.Name.DOS.Gaussian`、PDOS では `System.Name.PDOS.Tetrahedron.atom<I>.<軌道名>` などです。
  PDOS は原子ごと、\\( l \\) と多重度と \\( m \\) ごとの軌道（`s1`、`p1`（\\( p\_x \\)）など）と、その原子の和のファイルに分かれます。

### 3.4 `DosGauss.fileout`〔公開〕〔規約〕

〔公開〕マニュアル 23.2 節：固有ベクトルのファイルを書かずに、Gauss 型の広がりで DOS を直接求めます。四面体法は使えません。

〔規約〕`.Dos.val` は `mode 7` で、`Kgrid` は `1 1 1`、`<Eigenvalues` にはエネルギー格子の点が並び、`.Dos.vec` には格子点ごとに
Gauss 型で積算した軌道ごとの重みが入ります。DosMain は同じ手順で読み込みます。

## 4. 電荷解析

### 4.1 Mulliken 電荷〔公開〕〔規約〕

〔公開〕TN9 Eq. 1–6（記法を書き直した）：軌道 \\( a = (I\alpha) \\) の Mulliken の総和型の占有数（gross population）は

$$
M^{\sigma}\_{I\alpha} = \sum\_{\mathbf L}\sum\_{J\beta} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, S\_{I\alpha,J\beta}(\mathbf L)
$$

で、和則 \\( g\_s\sum\_\sigma\sum\_{I\alpha} M^\sigma\_{I\alpha} = N\_{\mathrm e} \\) を満たします（[密度行列と電荷ミキシング](density-mixing.md) の和則と同じ）。
原子の占有数は \\( M^\sigma\_I = \sum\_{\alpha\in I} M^\sigma\_{I\alpha} \\) です。k 空間で書けば、3.1 節の重み \\( w\_{a\mu} \\) を占有数で重み付けして足したものに等しくなります。

〔公開〕マニュアル 35.1 節：結果は既定で `System.Name.out` に書かれ、原子ごとの `Up spin`、`Down spin`、`Sum`、`Diff` と、
軌道ごとの分解（\\( m \\) についての和、\\( m \\) と多重度についての和を含む）が出ます。

〔規約〕

- 入力キーワードは不要で、SCF の各反復でも計算します（[非共線](noncollinear-energy-force.md) 6 節の原子のスピン、DFT+U の初期化などにも使う）。
- 出力の「電荷」は電子の占有数（正の数）です。原子の正味の電荷は \\( Z\_I^{\mathrm{val}} - (M^\uparrow\_I + M^\downarrow\_I) \\) で、これは cube ファイルの原子行（11.1 節）に書かれます。
- 非スピン偏極では `Up spin` と `Down spin` に同じ \\( M\_I \\)（1 スピン分）が出ます。
- 非共線では、2×2 の Mulliken 型の行列から \\( N^\uparrow, N^\downarrow, \theta, \varphi \\) を作ります（[非共線 DFT のエネルギーと力](noncollinear-energy-force.md) の 6.1 節）。
- 全スピン（`.scfout` の項目 20 などに使う）は共線で \\( S = \frac12\sum\_I (M^\uparrow\_I - M^\downarrow\_I) \\) です。
- 軌道の縮約を使った場合は、縮約後の基底での \\( \rho \\) と \\( S \\) を使います。

### 4.2 Voronoi 電荷（ファジーセル）〔公開〕〔規約〕

〔公開〕TN9 Eq. 7–16、マニュアル 35.2 節：Becke の fuzzy cell（Becke, JCP **88**, 2547 (1988)；Becke–Dickson, JCP **89**, 2993 (1988)）で空間を原子ごとに滑らかに分け、

$$
N^\sigma\_I = \int d^3r\\, w\_I(\mathbf r)\\, n^\sigma(\mathbf r),\qquad
w\_I(\mathbf r) = \frac{P\_I(\mathbf r)}{\sum\_J P\_J(\mathbf r)},\qquad
P\_I(\mathbf r) = \prod\_{J\ne I} s\left(\mu\_{IJ}(\mathbf r)\right)
$$

$$
\mu\_{IJ}(\mathbf r) = \frac{\lvert\mathbf r - \mathbf R\_I\rvert - \lvert\mathbf r - \mathbf R\_J\rvert}{\lvert\mathbf R\_I - \mathbf R\_J\rvert},\qquad
s(\mu) = \tfrac12\left[1 - p^{(k)}(\mu)\right],\qquad p(x) = \tfrac32 x - \tfrac12 x^3
$$

とします（\\( p^{(k)} \\) は \\( p \\) を \\( k \\) 回合成したもの）。原子の大きさによる補正はしないので、分割は構造だけで決まります。
キーワードは `Voronoi.charge`（`on`／`off`、既定 `off`）で、原子ごとの `Up spin`、`Down spin`、`Sum`、`Diff` とファジーセルの体積（Å³）が出ます。

〔規約〕

- **合成の回数は 6 回** です。TN9 は \\( k = 3 \\)（Becke の推奨値）と書いていますが、既存実装は \\( p \\) を合計 6 回合成します。
  回数が多いほど境界が鋭くなり、Voronoi 胞に近づきます。どちらの回数でも \\( \sum\_I w\_I = 1 \\) は成り立ちます（`tools/check_analysis.py`）。
- \\( P\_I \\) と分母の和の範囲：原子 \\( I \\) に属する格子点 \\( \mathbf r \\) について、分母の和は原子 \\( I \\) の近接リスト（1.2 節、自分自身を含む）の原子 \\( J \\)、
  各 \\( P\_J \\) の積は \\( J \\) の近接リスト（自分自身を除く）の原子についてとります。すべての原子についての和・積ではありません。
  近接リストは PAO のカットオフ半径で決まるので、厳密には「構造だけ」ではなくカットオフにも依存します。分母が \\( 10^{-14} \\) 未満なら \\( w\_I = 0 \\) とします。
- 積分は、原子 \\( I \\) の基底関数のカットオフ球に含まれる実空間格子点の和（格子点の体積を掛ける）です。\\( n^\sigma \\) は価電子密度（芯電荷を含まない）です。
  非共線では全密度と磁化ベクトルを \\( w\_I \\) で積分し、\\( N^\uparrow, N^\downarrow \\) と向きに直します。
- 体積は \\( \int w\_I\\,d^3r \\)（Å³ に換算）です。

### 4.3 静電ポテンシャルへのあてはめ（ESP 電荷）〔公開〕〔規約〕

〔公開〕TN9 Eq. 17–21、マニュアル 35.3 節：各原子に点電荷 \\( Q\_I \\) を置いた静電ポテンシャル
\\( V^{\mathrm{ESP}}(\mathbf r) = \sum\_I Q\_I/\lvert\mathbf r - \mathbf R\_I\rvert \\) を、DFT のポテンシャルに最小二乗であてはめます。
標本点は、原子を中心として van der Waals 半径の \\( \lambda\_1 \\) 倍と \\( \lambda\_2 \\) 倍の球の和集合で作る 2 つの殻の間の実空間格子点です。
全電荷の拘束（または全電荷と双極子の拘束）を Lagrange 乗数で課します（Singh–Kollman；Besler–Merz–Kollman）。
SCF で `level.of.fileout` を 1 以上にして `System.Name.out` と `System.Name.vhart.cube` を作り、プログラム esp を
`esp System.Name -c <0|1> -s λ1 λ2` として実行します。小分子向けで、大きな分子や結晶には向きません。

〔規約〕

- **あてはめる量**：`.vhart.cube` の値 \\( \delta V\_{\mathrm H}(\mathbf r) \\)（電子密度と原子密度の重ね合わせの差 \\( \delta n \\) の Hartree ポテンシャル、11.1 節）の **符号を反転したもの**
  \\( -\delta V\_{\mathrm H} \\) です。これは電荷 \\( -\delta n \\) の静電ポテンシャルで、中性原子の重ね合わせ（芯電荷と原子密度）の静電ポテンシャルは原子の外では 0 になるので、
  \\( Q\_I \\) は原子の正味の電荷（電子が増えれば負）になります。TN9 は「DFT の Hartree ポテンシャル」と書いていますが、正確にはこの差の量です。
- **方程式**：\\( b\_I(\mathbf r) = \sum\_{\mathbf n} 1/\lvert\mathbf r - \mathbf R\_I - \mathbf L\_{\mathbf n}\rvert \\)（周期像の和は cube ファイルの格子から決まる有限の範囲）、
  \\( A\_{IJ} = \sum\_p b\_I(\mathbf r\_p)\\,b\_J(\mathbf r\_p) \\)、\\( B\_I = -\sum\_p \delta V\_{\mathrm H}(\mathbf r\_p)\\,b\_I(\mathbf r\_p) \\) として、
  - `-c 0`：\\( \sum\_I Q\_I = 0 \\) を拘束（全電荷 \\( Q\_{\mathrm{tot}} \\) は **常に 0** で、`scf.system.charge` は使われない）。
  - `-c 1`：さらに \\( \sum\_I Q\_I \mathbf R\_I = \mathbf d \\)（\\( \mathbf d \\) は `.out` に書かれた全双極子を原子単位に直したもの）を拘束。

  の拡大した線形方程式を解きます。既定は `-c 0` です。
- **既定の倍率**：`-s` を省略すると \\( (\lambda\_1, \lambda\_2) = (1.4, 5.0) \\) です（マニュアルの例は `1.4 2.0`）。\\( \lambda\_2 > \lambda\_1 \ge 0 \\) でなければ誤りです。
- **van der Waals 半径**：Bondi の値（Å、例：H 1.20、C 1.70、N 1.55、O 1.52、F 1.47、Si 2.10、P 1.80、S 1.80、Cl 1.75）を原子番号で引きます。
  値の登録がない元素（遷移金属の多くなど）があると誤りで終了します。
- 原子の位置と格子は cube ファイルの見出しから読みます。双極子（Debye）は \\( \sum\_I Q\_I\mathbf R\_I \\) です。
- 出力の「RMS」は診断用の値で、定義が標準的な二乗平均平方根と異なる（各標本点・各原子・各周期像の差の二乗の和の平方根を、標本点数と格子点の体積で割る）ので、比較には使いません。
- 別の方式（`-l`：ファイルで与えた局所的な電荷に全電荷の拘束だけを課す）も残っていますが、マニュアルに記載がないので対象外とします。

## 5. 軌道分解全エネルギー〔公開〕〔規約〕

〔公開〕マニュアル 25 節：`Energy.Decomposition on`（既定 `off`）で、全エネルギーを原子と軌道に分解した結果を出します。

$$
U\_{\mathrm{tot}} = \sum\_\sigma\left(U\_{\mathrm{kin}}^\sigma + U\_{v}^\sigma + U\_{\mathrm{con}}^\sigma\right) + U\_{\mathrm{core}} + U\_{\mathrm{H0}} + U\_{\mathrm{vdW}}
$$

で、\\( U\_{\mathrm{con}} \\) は形式上加えた定数ポテンシャルの項です。実験的な機能とされています。

〔規約〕軌道 \\( a = (I\alpha) \\)、スピン \\( \sigma \\) ごとに

$$
U\_{\mathrm{kin},a}^\sigma = \sum\_{\mathbf L, b}\rho^\sigma\_{ab}(\mathbf L)\\,T\_{ab}(\mathbf L),\qquad
U\_{v,a}^\sigma = \sum\_{\mathbf L, b}\rho^\sigma\_{ab}(\mathbf L)\left[H^\sigma\_{ab}(\mathbf L) - T\_{ab}(\mathbf L)\right],\qquad
U\_{\mathrm{con},a}^\sigma = c\sum\_{\mathbf L, b}\rho^\sigma\_{ab}(\mathbf L)\\,S\_{ab}(\mathbf L)
$$

です（\\( T \\) は運動エネルギー行列、非スピン偏極では両スピン同じ）。定数 \\( c \\) は、全エネルギーから芯–芯、\\( E\_{\mathrm{H0}} \\)（中性原子の項）、
van der Waals の項を除いたものとバンドエネルギーの差を電子数で割った値で、和が全エネルギーに一致するように選ばれています
（各項の定義は [全エネルギーと力](energy-force.md)）。芯–芯と \\( E\_{\mathrm{H0}} \\) と vdW の項は原子ごとにだけ分け、その原子の最初の軌道に置きます。
非共線では \\( U\_v \\) に \\( \rho \\) と \\( H \\) の虚部・スピン非対角成分の積が加わります。
出力は `System.Name.DecE` で、原子ごとの表（\\( U\_{\mathrm{tot}} \\)、スピンごと、\\( U\_{\mathrm{kin}} \\)、\\( U\_v \\)、\\( U\_{\mathrm{con}} \\)、\\( U\_{\mathrm{core}} + U\_{\mathrm{H0}} \\)、\\( U\_{\mathrm{vdW}} \\)、Hartree）と軌道ごとの表です。
軌道の縮約とは同時に使えません。

## 6. 交換結合定数（jx）

### 6.1 公式〔公開〕

Liechtenstein らの磁気力の定理（JMMM **67**, 65 (1987)）に基づき、Han–Ozaki–Yu（PRB **70**, 184421 (2004)）が非直交の局在基底に拡張した公式です。
原子 \\( i \\) の交換分裂の行列を

$$
\Delta\_i = H^{\uparrow}\_{ii}(\mathbf 0) - H^{\downarrow}\_{ii}(\mathbf 0)\qquad(n\_{\mathrm{orb}}(i)\times n\_{\mathrm{orb}}(i)\text{ の一中心ブロック})
$$

とし、固有ベクトル \\( c^\sigma\_{\mu}(\mathbf k) \\) の原子 \\( i \\)、\\( j \\) の部分を使って

$$
J\_{ij} = \frac14\\,\frac{1}{N\_k}\sum\_{\mathbf k}\sum\_{\mu\nu}
\frac{f(\varepsilon^\uparrow\_\mu) - f(\varepsilon^\downarrow\_\nu)}{\varepsilon^\downarrow\_\nu - \varepsilon^\uparrow\_\mu}\\,
\left[c^{\uparrow\dagger}\_{\mu}\Delta\_j c^{\downarrow}\_{\nu}\right]\left[c^{\downarrow\dagger}\_{\nu}\Delta\_i c^{\uparrow}\_{\mu}\right]
$$

とします（マニュアル 45.1 節 Eq. 42、47 に対応。\\( \mathbf k \\) の引数は省略）。これは周期像について和をとった \\( \sum\_{\mathbf L} J\_{i\mathbf 0, j\mathbf L} \\) です。
個々の \\( J\_{i\mathbf 0,j\mathbf L} \\) は Green 関数

$$
G^\sigma\_{ij}(z;\mathbf L) = \frac{1}{N\_k}\sum\_{\mathbf k} e^{-i\mathbf k\cdot\mathbf L}\sum\_\mu \frac{c^\sigma\_{i\mu}(\mathbf k)\\,c^{\sigma\dagger}\_{j\mu}(\mathbf k)}{z - \varepsilon^\sigma\_\mu(\mathbf k)}
\qquad\left(= \left[(zS - H^\sigma)^{-1}\right]\_{i\mathbf 0, j\mathbf L}\right)
$$

と Fermi 関数の連分数展開の極（Ozaki, PRB **75**, 035123 (2007)；Terasawa ら, JPSJ **88**, 114706 (2019)）を使って

$$
J\_{i\mathbf 0, j\mathbf L} = \frac{k\_{\mathrm B}T}{2}\sum\_{p=1}^{N\_P} R\_p\\, \mathrm{Re}\\,\mathrm{Tr}\left[\Delta\_i\\, G^\downarrow\_{ij}(z\_p;\mathbf L)\\, \Delta\_j\\, G^\uparrow\_{ji}(z\_p;-\mathbf L)\right],\qquad
z\_p = \mu + i\\,k\_{\mathrm B}T\\,\zeta\_p
$$

で求めます（マニュアル Eq. 44–46）。極 \\( \zeta\_p \\) と留数 \\( R\_p \\) は、\\( M = 2N\_P \\) 次の一般化固有値問題 \\( A\mathbf x = \lambda B\mathbf x \\)
（\\( A \\) は対角が 0、隣接成分が \\( -\tfrac12 \\) の三重対角、\\( B = \mathrm{diag}(1, 3, 5, \dots, 2M-1) \\)、\\( \mathbf x^{\mathsf T}B\mathbf x = 1 \\)）の
負の固有値 \\( \lambda\_p \\)（小さい方から \\( N\_P \\) 個）と固有ベクトルの第 1 成分 \\( x\_{p,1} \\) から

$$
\zeta\_p = -1/\lambda\_p,\qquad R\_p = -\tfrac14\\, x\_{p,1}^2\\, \zeta\_p^2
$$

で、\\( 1/(1+e^{x}) = \tfrac12 + \sum\_p R\_p\left[(x - i\zeta\_p)^{-1} + (x + i\zeta\_p)^{-1}\right] \\) を与えます
（`tools/check_analysis.py`：\\( N\_P = 60 \\) で \\( \lvert x\rvert \le 30 \\) の誤差 \\( 10^{-15} \\)。低い極は松原の極 \\( \zeta = (2p-1)\pi \\)、\\( R\_p = -1 \\) に一致）。
\\( \sum\_{\mathbf L} J\_{i\mathbf 0,j\mathbf L} \\) は上の周期和の公式に一致します（同スクリプトで \\( 10^{-17} \\)）。

### 6.2 符号と規格化の規約〔設計〕（検算）

磁気力の定理（固定した \\( \mu \\) での大正準ポテンシャル \\( \Omega \\) の、スピンの向きの回転に対する 2 次の変化）と比べると、
上の \\( \frac14 \\) の公式（Band のモード）の \\( J \\) は

$$
{}-\frac{\partial^2\Omega}{\partial\theta\_i\\,\partial\theta\_j}\bigg\rvert\_{\theta=0} = 2J\_{ij}
\quad\Longleftrightarrow\quad
E = -\sum\_{i\ne j} J\_{ij}\\,\mathbf e\_i\cdot\mathbf e\_j = -2\sum\_{i<j}J\_{ij}\\,\mathbf e\_i\cdot\mathbf e\_j
$$

に対応します（\\( \mathbf e\_i \\) はスピンの向きの単位ベクトル。`tools/check_analysis.py` の 3 サイトの非直交クラスターで \\( 10^{-8} \\) の精度で確認）。
\\( J > 0 \\) が強磁性的です。スピンの大きさは \\( J \\) に繰り込まれています。

### 6.3 プログラムと設定ファイル〔公開〕〔規約〕

〔公開〕マニュアル 45.2–45.6 節：SCF を `HS.fileout on` で行い、`jx System.Name.scfout jx.config` と実行します。
共線計算（Cluster または Band）だけに対応し、非共線は対象外です。

〔規約〕設定ファイルはキーワード形式で、次の項目を読みます。

| キーワード | 既定値 | 内容 |
|---|---|---|
| `Flag.PeriodicSum` | `off` | Band のとき、`on` なら周期和 \\( \sum\_{\mathbf L}J\_{i\mathbf 0,j\mathbf L} \\)、`off` なら個々の \\( J\_{i\mathbf 0,j\mathbf L} \\) |
| `Num.Poles` | 60 | 極の数 \\( N\_P \\)（`PeriodicSum off` でのみ使う） |
| `Num.Kgrid` | `1 1 1` | k 点の分割数（0 以下は誤り） |
| `Num.ij.pairs` | 0 | 計算する対の数 |
| `Bunch.ij.pairs` | `Num.ij.pairs` | 一度にまとめて計算する対の数（結果に影響しない） |
| `<ijpairs.cellid` … | — | 各行 `i j n1 n2 n3`（原子番号は 1 始まり、\\( \mathbf L = \sum n\_l\mathbf a\_l \\)） |
| `<ijpairs.nocellid` … | — | 各行 `i j` |
| `Flag.LNO` | `off` | 局在自然軌道（LNO）で \\( \Delta \\) を射影する（Band のみ） |
| `LNO.occ.cutoff` | 0.0 | LNO の占有数の閾値 |

- Band で `PeriodicSum off` のときは `<ijpairs.cellid` だけを読みます。それ以外は `<ijpairs.nocellid` を優先し、なければ `<ijpairs.cellid` の格子ベクトルを無視して使います。
- **解法の振り分け**：`.scfout` の解法の番号が Cluster（2）または Cluster-DIIS（7）なら Cluster のモード（Γ 点、MPI 1 プロセスのみ）、Band（3）なら Band のモード、
  それ以外では何も計算しません。
- **k 点**：\\( q\_{i,l} = -\tfrac12 + (2l+1)/(2N\_i) \\)（Monkhorst–Pack 型、\\( N\_i \\) が偶数なら Γ を含まない）、すべて重み \\( 1/N\_k \\)、削減なし。
- **固有値問題**：\\( S(\mathbf k) \\) の固有分解による直交化で、\\( S \\) の小さな固有値の正則化はしません（[一般化固有値問題](generalized-eigen.md) のどの規則とも異なる）。
- **温度と化学ポテンシャル**：`.scfout` の \\( \mu \\) と \\( T \\)、\\( k\_{\mathrm B} = 3.166813628\times10^{-6} \\) Hartree/K です。
- **個々の \\( J \\) の数値**：\\( z\_p \\) の虚部に \\( \eta = 10^{-10} \\) Hartree を加えます（影響は無視できる）。
- **Cluster のモード**：Γ 点だけで同じ形の和をとりますが、**値は Band のモードの公式の \\( \tfrac12 \\) 倍** です
  （交換分裂を \\( \tfrac12\Delta \\) とし、前の係数を \\( \tfrac12 \\) とするため、全体で \\( \tfrac18 \\)）。同じ Γ 点の問題を Band のモード（`Num.Kgrid 1 1 1`、周期和）で解くと 2 倍の値になります。
  どちらが意図した規格化かは公開資料から判断できません（〔未確定〕）。新しい実装は Band のモードの規格化（6.2 節）を既定とし、Cluster の規格化を互換オプションとします。
- **LNO**：原子ごと・スピンごとに \\( (\rho S)\_{ii} \\) の一中心ブロック（相手の原子について和をとったもの）を対角化し、占有数が閾値を超える固有ベクトルへの射影 \\( P^\sigma\_i \\) で
  \\( \Delta\_i \to P^\uparrow\_i\Delta\_iP^\downarrow\_i \\)（\\( j \\) 側は \\( P^\downarrow\_j\Delta\_jP^\uparrow\_j \\)）と置き換えます。固有ベクトルは非対称行列の右固有ベクトルです。
- **出力**：標準出力だけで、各対について \\( J \\) を meV と mRy（\\( 1\\,\mathrm{Ha} = 2000\\,\mathrm{mRy} \\)）で出します。Band で `PeriodicSum off` のときは格子ベクトルの整数座標も出します。
- `i = j` の周期和には \\( \mathbf L = \mathbf 0 \\)（自分自身）の項が含まれます。物理的な結合として使うときは個々の \\( J \\) の和から除きます。

## 7. Berry 位相による電気分極（polB）

### 7.1 公式〔公開〕

King-Smith–Vanderbilt（PRB **47**, 1651 (1993)）と Resta（RMP **66**, 899 (1994)）の Berry 位相の公式を、TN11 Eq. 1–11 に従って LCAO で書きます。
逆格子ベクトル \\( \mathbf b\_1 \\) の方向について、\\( \mathbf b\_1 \\) に沿った k 点の列（弦）\\( \mathbf k\_j = \mathbf k\_\perp + (j/N\_1)\mathbf b\_1 \\)、\\( j = 0, \dots, N\_1 \\) と
\\( \mathbf b = \mathbf b\_1/N\_1 \\) について、占有状態の間の重なり行列

$$
M\_{\mu\nu}(\mathbf k, \mathbf k + \mathbf b) = \langle\psi\_{\mu\mathbf k}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\psi\_{\nu,\mathbf k+\mathbf b}\rangle\quad(\text{単位胞あたり})
$$

の行列式の積の偏角から、弦ごとの Berry 位相

$$
\varphi(\mathbf k\_\perp) = -\arg\prod\_{j=0}^{N\_1-1}\det M(\mathbf k\_j, \mathbf k\_{j+1})
$$

を求め、弦について平均した \\( \bar\varphi\_1 \\) から、電子の寄与の双極子（単位胞あたり、\\( e > 0 \\)、電子の電荷 \\( -e \\)）

$$
\mathbf d\_{\mathrm{el}} = -\frac{g\_s\\,e}{2\pi}\sum\_{i=1}^{3}\bar\varphi\_i\\,\mathbf a\_i,\qquad \mathbf P = \frac{\mathbf d\_{\mathrm{el}} + \mathbf d\_{\mathrm{ion}}}{\Omega}
$$

を得ます。共線スピン偏極ではスピンごとの \\( \varphi \\) を足し、\\( g\_s = 1 \\) とします。

LCAO では、Bloch 和が格子ベクトルだけの位相なので \\( c(\mathbf k + \mathbf G) = c(\mathbf k) \\) と選べ、弦の端 \\( \mathbf k\_{N\_1} = \mathbf k\_0 + \mathbf b\_1 \\) の係数として \\( \mathbf k\_0 \\) の係数をそのまま使えます
（周期的ゲージ）。行列要素は

$$
M\_{\mu\nu}(\mathbf k, \mathbf k + \mathbf b) = \sum\_{ab}\overline{c\_{a\mu}(\mathbf k)}\\,T\_{ab}\\,c\_{b\nu}(\mathbf k + \mathbf b),\qquad
T\_{I\alpha,J\beta} = \sum\_{\mathbf L} e^{i(\mathbf k+\mathbf b)\cdot\mathbf L}\\,\langle\phi\_{I\alpha}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\phi\_{J\beta}(\mathbf L)\rangle
$$

で、TN11 Eq. 10–11 は \\( e^{-i\mathbf b\cdot\mathbf r} \\) を行側の原子の中心のまわりで 1 次まで展開します。

$$
\langle\phi\_{I\alpha}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\phi\_{J\beta}(\mathbf L)\rangle \approx X\_{I\alpha,J\beta}(\mathbf L) \equiv
e^{-i\mathbf b\cdot\mathbf R\_I}\left[S\_{I\alpha,J\beta}(\mathbf L) - i\\,\mathbf b\cdot\mathbf D\_{I\alpha,J\beta}(\mathbf L)\right]
$$

（\\( \mathbf D \\) は `.scfout` の項目 15 の位置演算子の行列）。誤差は \\( O(\lvert\mathbf b\rvert^2) \\) で、弦の分割を細かくすると消えます。
`tools/check_analysis.py` は、点状の軌道の Rice–Mele 鎖でこの構成（格子ベクトルだけの位相と周期的ゲージ）が Resta の超格子の公式と一致することを確認しています。

### 7.2 既存実装の選択〔規約〕

- **入力**：`polB System.Name.scfout` を実行し、k 点の分割数 \\( N\_1\\; N\_2\\; N\_3 \\)（各 1 以上）と、計算する方向（逆格子の 3 方向それぞれ 0 か 1）を順に入力します（マニュアル 44 節）。
  解法とスピンの種類（0、1、3）の制限はありません。
- **弦**：方向 \\( i \\) について、他の 2 方向の分数座標 \\( l/N \\)（\\( l = 0, \dots, N - 1 \\)、Γ 中心）のすべての組で弦を作り、弦の上の点も \\( j/N\_i \\) です。
  弦の数は他の 2 方向の分割数の積です。
- **\\( T \\) の対称化**：上の \\( X \\) から、\\( T \\) を 2 通りに作って平均します（和はすべての近接対 \\( (I, J, \mathbf L) \\)）。

  $$
  T = \tfrac12\left(T^{(1)} + T^{(2)}\right),\qquad
  T^{(1)}\_{I\alpha,J\beta} = \sum\_{\mathbf L} e^{i(\mathbf k+\mathbf b)\cdot\mathbf L}X\_{I\alpha,J\beta}(\mathbf L),\qquad
  T^{(2)}\_{J\beta,I\alpha} = \sum\_{\mathbf L} e^{-i\mathbf k\cdot\mathbf L}X\_{I\alpha,J\beta}(\mathbf L)
  $$

  \\( T^{(2)} \\) は \\( \langle\phi\_{J\beta}\rvert e^{-i\mathbf b\cdot\mathbf r}\lvert\phi\_{I\alpha}(-\mathbf L)\rangle \\) を、\\( J \\) ではなく \\( I \\) の中心のまわりの展開で近似したもので
  （\\( S \\) と \\( \mathbf D \\) が実数なので複素共役をとらない形になる）、2 つの平均は 1 次の展開の中心の選び方の偏りを打ち消します。
- **占有状態の数**：非スピン偏極で \\( \lfloor N\_{\mathrm e}/2\rfloor \\)、非共線で \\( \lfloor N\_{\mathrm e}\rfloor \\)、共線スピン偏極ではスピンごとに、その k 点で \\( \varepsilon \le \mu \\) の状態の数です。
  非スピン偏極・非共線で、この数が \\( \varepsilon \le \mu \\) の状態の数と異なる k 点があれば「金属」と警告します（計算は続ける）。
- **固有値問題**：一般化固有値問題をそのまま解きます（正則化なし）。非共線では \\( 2N \\) 次のスピノル行列で、\\( T \\) は両スピンに同じものを使います。
- **弦の平均**：各弦の \\( \psi = \arg\prod\det M \in (-\pi, \pi] \\) について、基準 \\( \psi\_0 = \arg\sum\_{\text{弦}} e^{i\psi} \\) を求め、
  各 \\( \psi \\) を \\( 2\pi \\) の整数倍だけずらして \\( (\psi\_0 - \pi, \psi\_0 + \pi] \\) に入れてから平均します（\\( \bar\varphi = -\overline{\psi} \\)）。
- **出力**：電子・芯（イオン）・背景電荷の寄与と合計の双極子（Debye）と分極（μC/cm²）を、Cartesian の 3 成分で標準出力に出します。
  芯と背景電荷の寄与は `.scfout` の項目 20 の値をそのまま使います（芯の寄与を \\( 2\pi \\) の任意性で折り返すことはしない）。
  換算は 1 \\( e \\)·bohr = 2.54174776 Debye、\\( P[\mu\mathrm C/\mathrm{cm}^2] = 5721.52891433\\,d[e\\,\mathrm{bohr}]/\Omega[\mathrm{bohr}^3] \\) です。
  分極の量子として \\( e\mathbf a\_i/\Omega \\) も出しますが、非スピン偏極では実際の任意性は \\( g\_s e\mathbf a\_i/\Omega \\)（2 倍）です。
- 分極は分極の量子の整数倍の任意性を持つので、意味があるのは差（Born 有効電荷 \\( Z^\*\_{\alpha\beta} = (\Omega/e)\\,\Delta P\_\alpha/\Delta u\_\beta \\) など）です（マニュアル 44 節）。

### 7.3 SCF 中の双極子〔規約〕

孤立系の双極子 \\( \sum\_I Z\_I^{\mathrm{val}}\mathbf R\_I - \int n(\mathbf r)\\,\mathbf r\\,d^3r \\)（背景電荷の項を含む）は、SCF の後に格子上の和で求めて `.out` に Debye で出します。
既存実装には SCF の中で Berry 接続の解析式を使う分極の経路もありますが、入力キーワードとの対応がマニュアルに記載されていないので〔未確定〕とします。

## 8. Chern 数・Berry 曲率・\\( \mathbb Z\_2 \\) 不変量

### 8.1 Fukui–Hatsugai–Suzuki の方法〔公開〕

Fukui–Hatsugai–Suzuki（JPSJ **74**, 1674 (2005)）とマニュアル 58 節：2 次元の k 面を格子に分け、隣り合う k 点の間のリンク変数
\\( U(\mathbf k, \mathbf k') = \det M(\mathbf k, \mathbf k')/\lvert\det M\rvert \\)（\\( M \\) は 7.1 節と同じ重なり行列、選んだバンドについて）から、
小区画（プラケット）ごとの格子場の強さ

$$
F = \arg\left[U(\mathbf k\_1, \mathbf k\_2)\\,U(\mathbf k\_2, \mathbf k\_3)\\,U(\mathbf k\_3, \mathbf k\_4)\\,U(\mathbf k\_4, \mathbf k\_1)\right]\in(-\pi, \pi]
$$

を求め、Chern 数 \\( C = \frac{1}{2\pi}\sum F \\) を得ます。メッシュが十分に細かければ \\( C \\) は整数で、Berry 曲率は \\( F \\) を小区画の面積で割ったものです。
\\( \mathbb Z\_2 \\) 不変量（Fukui–Hatsugai, JPSJ **76**, 053702 (2007)；Sawahata ら, JJAP **57**, 030309 (2018)；マニュアル 59 節）は、
時間反転不変な面の半分の Brillouin 帯で、時間反転不変運動量（TRIM）で Kramers 対のゲージを固定したうえで、
\\( D = \frac{1}{2\pi}\left[\sum\_{\text{半分}} F - \oint A\right] \bmod 2 \\) として求めます。

### 8.2 既存実装（calB、Z2FH）〔規約〕

- どちらも非共線（\\( s = 3 \\)）の `.scfout` だけを受け付けます。重なり行列は 7.2 節と同じ 1 次展開と対称化で作ります。
- **calB**：対話的に (a) 1 バンドの Berry 曲率（\\( M \\) の対角成分）か、指定した数までのバンドの和（行列式）か、(b) バンドの数（0 以下なら価電子数）、
  (c) 面の法線の方向（逆格子の 3 方向の 0／1）、(d) 面内の分割数 2 つ、を読みます。面は Γ を通り、\\( [0,1)\times[0,1) \\) を等分します。
  出力 `BerryC<バンド数>.dat` に分割数、バンド数、Chern 数の見出しと、各小区画の中心の分数座標と、\\( F \\) を面積（Å²）で割った値を書きます。
  面積として、面内の 2 つの逆格子ベクトルの外積の、法線方向の **Cartesian 成分**（大きさではない）を使います。また 2 列目の座標は、2 つの分割数が異なると誤った値になります。
  Chern 数そのものはこれらの影響を受けません。
- **Z2FH**：対話的に (a) 半分の方向の分割数 \\( m \\)、(b) \\( k = 0 \\) の面を計算するかどうか、(c) 開始する面（6 面：各方向の \\( 0 \\) と \\( \pi \\)）を読みます。
  面の中の 2 方向を \\( [-\tfrac12, \tfrac12] \\)（\\( 2m \\) 分割）と \\( [0, \tfrac12] \\)（\\( m \\) 分割）に分け、TRIM の 4 点で
  \\( \Theta(\psi\_\uparrow, \psi\_\downarrow) = (\overline{\psi\_\downarrow}, -\overline{\psi\_\uparrow}) \\) により Kramers 対の一方から他方を作り、半分の境界上の状態も \\( \Theta \\) で作ります。
  小区画ごとの \\( n = \frac{1}{2\pi}\left[\sum\_{\text{辺}}\arg U - \arg\prod U\right] \\) の和を mod 2 にし、6 面の値から \\( (\nu\_0; \nu\_1\nu\_2\nu\_3) \\) を `Z2.dat` に出します。
  面内の第 2 方向（分割 \\( m \\)、幅 \\( \tfrac12 \\)）のリンクでは、\\( e^{-i\mathbf b\cdot\mathbf r} \\) の展開に使う \\( \mathbf b \\) が実際の刻み \\( 1/(2m) \\) ではなく \\( 1/m \\) です。
  Bloch 位相は正しい k を使うので、影響は位置演算子の補正の大きさだけです。整数の結果が変わる条件は〔未確定〕で、新しい実装は正しい \\( \mathbf b \\) を使います。

## 9. k 空間のスピンテクスチャ（kSpin）〔公開〕〔規約〕

〔公開〕マニュアル 56 節、技術ノート「Analysis of spin texture in the k-space」、Kotaka ら（JJAP **52**, 035204 (2013)）：
非共線の Band 計算の `.scfout` から、指定した k 点の各状態について、原子・軌道ごとに分解したスピンの期待値を求めます。
等エネルギー線の追跡（FermiLoop）、面上の格子（GridCalc）、経路に沿ったバンド（BandDispersion）、指定した k 点（MulPOnly）のモードと、
結果を原子について集計する別プログラム（MulPCalc）があります。

〔規約〕

- 入力はキーワード形式（`Filename.scfout`、`Filename.outdata`、`Calc.Type`、`Energy.Range`（eV、\\( \mu \\) 基準、既定 \\( -0.5\dots0.5 \\)）、
  `Band.Nkpath` と `<Band.kpath`、`Search.kCentral`、`Calc.Type.3mesh`、`k-plane.1stStep`、`k-plane.2ndStep`、`kRange.3mesh`、`k-plane.height`、`kRange.height`、
  `MulP.Vec.Scale`、`Spin.Degenerate`、`Calc.Bandbyband` など）です。解法は Band、スピンは非共線に限ります。
- 状態 \\( \mu \\) の軌道 \\( a \\) ごとのスピン密度行列を、Mulliken の総和型

  $$
  n^{\sigma\sigma'}\_{a\mu} = \overline{c\_{\sigma a,\mu}}\left(S(\mathbf k)\\,c\_{\sigma'\mu}\right)\_a
  $$

  で作り、原子（または軌道）について和をとって \\( \mathbf m = (2\\,\mathrm{Re}\\,n^{\uparrow\downarrow}, 2\\,\mathrm{Im}\\,n^{\uparrow\downarrow}, n^{\uparrow\uparrow} - n^{\downarrow\downarrow}) \\)
  （Pauli 行列の期待値の単位、\\( \hbar/2 \\) 倍ではない）とし、`MulP.Vec.Scale` を掛けて出力します。全軌道の和は \\( \langle\psi\rvert\boldsymbol\sigma\lvert\psi\rangle \\) に一致します。
- エネルギーは \\( \varepsilon - \mu \\) を eV で出し、k 点は Cartesian（bohr⁻¹、MulPCalc では Å⁻¹）で出します。
- FermiLoop は粗い格子で符号の変わる小区画を探し、細かい格子と Brent 法（既定で有効、反復 5 回）で \\( \varepsilon(\mathbf k) = E \\)（\\( E \\) は `Energy.Range` の中点）の点を求めて線をつなぎます。

## 10. 光学伝導度と誘電関数（CDDF）〔公開〕〔規約〕

〔公開〕マニュアル 62 節、技術ノート「Calculation of Optical Matrix Elements」「Dielectric Function and Conductivity」：
SCF の後に Kubo–Greenwood 型の式で光学伝導度テンソル \\( \sigma\_{ab}(\omega) \\) と誘電関数テンソル \\( \varepsilon\_{ab}(\omega) \\) を求めます。
運動量の行列要素は LCAO の \\( \langle\phi\rvert\nabla\lvert\phi'\rangle \\)（方法 0）、または \\( [H, \mathbf r] \\) による形（方法 1、非局所擬ポテンシャルの寄与を含む、非共線にも対応）から作ります。

〔規約〕キーワードと既定値：

| キーワード | 既定値 | 内容 |
|---|---|---|
| `CDDF.start` | `off` | |
| `CDDF.FWHM` | 0.2 eV | 広がり \\( \eta \\)（Lorentz 型の分母の \\( i\eta \\) に **FWHM の値そのもの** を使う） |
| `CDDF.maximum_energy` | 10.0 eV | 出力する周波数の上限 |
| `CDDF.minimum_energy` | 0.0 eV | 刻み \\( (E\_{\max} - E\_{\min})/N \\) にだけ影響（格子の始点は動かない） |
| `CDDF.additional_maximum_energy` | 0.0 eV | 含めるバンドの範囲を広げる |
| `CDDF.frequency.grid.total_number` | 10000 | \\( N \\) |
| `CDDF.approach` | 0 | 0：\\( \langle\phi\rvert\nabla\lvert\phi'\rangle \\)（共線のみ）、1：\\( [H,\mathbf r] \\) の形（共線・非共線） |
| `CDDF.shift_ChemP` | 0.0 eV | \\( \mu \\) をずらす |
| `CDDF.Kgrid` | `scf.Kgrid` の値 | `scf.Kgrid` が `1 1 1` のときに限り `1 1 1` でなければならない |

- 周波数は \\( \omega\_j = j\\,\Delta \\)（\\( j \ge 1 \\)）で、\\( \omega = 0 \\) の値は \\( j = 1, 2 \\) からの線形外挿 \\( 2y\_1 - y\_2 \\) です。
- k 点は Monkhorst–Pack 型 \\( -\tfrac12 + (2l+1)/(2N) \\)（\\( N = 1 \\) なら 0）で、\\( \pm\mathbf k \\) を重み 2 の 1 点にまとめます。
- 方法 0 の運動量の行列要素は、Bloch 和を格子ベクトルだけの位相でとった \\( \nabla^a(\mathbf k) \\) から \\( p^{a}\_{mn}(\mathbf k) \propto c^\dagger\_m\\,\nabla^a(\mathbf k)\\,c\_n \\) です。
- 帯間の項は、\\( \Delta\varepsilon = \varepsilon\_m - \varepsilon\_n \\)、\\( \Delta f = f\_n - f\_m \\)（非スピン偏極では占有数を 2 倍）として
  \\( \sigma\_{ab} \propto \sum \frac{\Delta f}{\Delta\varepsilon}\\, p^a\_{mn}\overline{p^b\_{mn}}\\,\frac{i}{\Delta\varepsilon + \omega + i\eta} \\)、
  \\( \varepsilon\_{ab} - \delta\_{ab} \propto \sum \frac{-i}{\Delta\varepsilon}\times(\text{同じ項}) \\) です。係数は \\( e^2/(\hbar a\_0) = 4599847.99 \\) S/m を \\( \Omega\sum w\_{\mathbf k} \\) で割ったもの（\\( \sigma \\)）、
  さらに \\( \varepsilon\_0 E\_{\mathrm h}/\hbar = 366044.288 \\) S/m で割ったもの（\\( \varepsilon \\)）です。
  帯内（Drude）の項は、非共線の経路では \\( -f'(\varepsilon) \\) を使う形で入りますが、共線の経路では \\( \Delta f = 0 \\) の対を最初に除くので入りません。
- 出力：`System.Name.cd_re`、`.cd_im`（\\( \sigma \\)、S/m）、`.df_re`、`.df_im`（\\( \varepsilon \\)。**\\( \delta\_{ab} \\) を含まない**）。
  1 列目がエネルギー（eV）、2–10 列目が \\( xx, xy, xz, yx, yy, yz, zx, zy, zz \\)、11 列目が対角和の 1/3、共線スピン偏極ではスピンごとの 9 成分が続きます。
  派生量（屈折率、消衰係数、吸収係数、反射率、透過率）のファイルでは \\( \delta\_{ab} \\) を加えます。
- \\( \sigma \\) の位相（\\( i \\) の向き）と \\( p^a\overline{p^b} \\) の添字の順は非対角成分の符号に影響します。公開資料の式との照合は〔未確定〕です。

## 11. 格子上の量の出力

### 11.1 cube ファイル〔公開〕〔規約〕

〔公開〕マニュアル 9・21 節：`level.of.fileout`（0–3、既定 1）が 1 以上なら、電子密度とポテンシャルを Gaussian cube 形式で出力します。

〔規約〕

| ファイル | 内容 | 条件 |
|---|---|---|
| `.tden.cube` | 全電子密度 \\( n^\uparrow + n^\downarrow \\)（価電子、bohr⁻³） | level ≥ 1 |
| `.den0.cube`、`.den1.cube` | 上向き・下向きの密度 | スピン偏極・非共線 |
| `.sden.cube` | \\( n^\uparrow - n^\downarrow \\) | 同上 |
| `.dden.cube` | \\( n - n\_{\mathrm{atom}} \\)（原子密度の重ね合わせとの差） | level ≥ 1 |
| `.v0.cube`、`.v1.cube` | スピンごとの有効ポテンシャル \\( \delta V\_{\mathrm H} + V\_{\mathrm{xc},\sigma} \\)（＋外部電場）（Hartree）。`scf.ProExpn.VNA` が `off` のときだけ \\( V\_{\mathrm{NA}} \\) を含む（既定の `on` では含まない） | level ≥ 1 |
| `.vhart.cube` | \\( \delta V\_{\mathrm H} \\)：\\( \delta n = n - n\_{\mathrm{atom}} \\) の Hartree ポテンシャル（\\( \mathbf G = 0 \\) の成分は 0）（Hartree） | level ≥ 1 |
| `.vna.cube`、`.vxc0.cube`、`.vxc1.cube`、`.grid` | 中性原子ポテンシャル、交換相関ポテンシャル、格子点の座標（Å） | level ≥ 2 |
| `.ncsden.xsf`、`.nc.xsf`、`.nco.xsf` | 非共線のスピン密度ベクトル場、原子ごとのスピン・軌道モーメント（XCrySDen） | 非共線 |

- **見出し**：1–2 行目は説明（分子軌道では固有値と \\( \mu \\)）、3 行目は原子数と格子の原点、4–6 行目は各方向の格子点数と格子の刻みのベクトル（bohr、点数は正）、
  続いて原子ごとに原子番号、正味の電荷（\\( Z^{\mathrm{val}} \\) から Mulliken の電子数を引いた値）、位置（bohr）です。
- **データ**：第 1 方向が最も外側、第 3 方向が最も内側（最も速く変わる）で、1 行に 6 個、第 3 方向の 1 列ごとに改行します。
- `.v0` の \\( V\_{\mathrm{NA}} \\) は短距離なので真空中では 0 です。ポテンシャルの基準（定数の任意性）は \\( \delta V\_{\mathrm H} \\) の \\( \mathbf G = 0 \\) の扱いで決まり、
  ESM（[ESM・電場・内殻正孔](esm-field-corehole.md)）を使う場合は真空の値が意味を持ちます。

### 11.2 分子軌道と LCAO 係数〔公開〕〔規約〕

〔公開〕マニュアル 33・34 節：Cluster では `level.of.fileout 2` で `.out` に LCAO 係数と固有値（Hartree）が出ます。
`MO.fileout on` と `num.HOMOs`、`num.LUMOs`、周期系では `MO.Nkpoint` と `<MO.kpoint`（分数座標）で、HOMO・LUMO の付近の分子軌道を cube ファイルに出力します。

〔規約〕`num.HOMOs`、`num.LUMOs`、`MO.Nkpoint` の既定値はいずれも 1 です。Cluster、Band、Cluster-DIIS でのみ有効です。HOMO は占有数が \\( \tfrac12 \\) を超える最後の状態です。
\\( \psi(\mathbf r) = \sum\_a c\_a\phi\_a(\mathbf r) \\)（\\( c^\dagger Sc = 1 \\)、追加の規格化なし、単位 bohr^{-3/2}）の実部・虚部を
`System.Name.homo<k>_<spin>_<番号>_r.cube`、`_i.cube`（Cluster では k の番号なし。虚部はスピン軌道相互作用や非共線のときだけ）、LUMO も同様に書きます。
周期系では、格子点ごとに基底関数の周期像 \\( \mathbf L \\) に Bloch 位相を付けて足します。

### 11.3 STM 像のための部分電荷密度〔公開〕〔規約〕

〔公開〕マニュアル 53 節：Tersoff–Hamann の近似で、\\( \mu \\) から測ったエネルギーの窓の中の状態の電荷密度を `System.Name.pden.cube` に出します。
`partial.charge on`、`partial.charge.energy.window`（eV、正なら伝導帯側、負なら価電子帯側）と、`Dos.fileout`（または `DosGauss.fileout`）が必要です。

〔規約〕窓の既定値は 0.0 eV です。Cluster と Band でのみ使えます。状態ごとの重み \\( \lvert f(\varepsilon) - f(\varepsilon - \Delta E)\rvert \\)（\\( \Delta E \\) は窓、\\( f \\) は計算の温度での Fermi 分布）で
密度行列を作り、格子上の密度に直します。非スピン偏極では 2 倍し、スピン偏極では `.pden0.cube`、`.pden1.cube` も出します。

### 11.4 仕事関数〔公開〕〔規約〕

〔公開〕マニュアル 50 節：スラブ模型で、面内で平均した静電ポテンシャルの真空での値 \\( \phi\_\infty \\) と \\( \mu \\) から \\( W = \phi\_\infty - \mu \\) とします
（\\( E[N-1] - E[N] = -\mu \\) による。ギャップのある系では近似）。分極した表面では ESM が必要です。

〔規約〕面内平均には補助プログラム gcube2oned を使い、cube ファイル（通常は `.v0.cube`）と軸の番号を与えると、
格子点の番号、軸に沿った位置（Å）、面内の平均値（cube の単位のまま）を出します。
軸 2・3 を選んだとき、出力の 1 列目（番号）は正しく進みません（2 列目の位置は正しい）。

## 実装への要求〔設計〕

- **入力は疎ブロック形式**：すべての解析は、構造（原子位置、格子、近接リスト）と、[疎ブロック形式](../design/sparse-block-format.md) の \\( H^\sigma(\mathbf L) \\)、\\( S(\mathbf L) \\)、
  \\( \rho^\sigma(\mathbf L) \\)、必要なら \\( \mathbf D(\mathbf L) \\)、\\( \nabla(\mathbf L) \\) を入力とする、状態を持たない関数にします。`.scfout` はそれらへの変換器（読み手・書き手）として 1 か所に閉じ込めます。
- **共通部品**：Bloch 和、k 点の生成（Γ 中心・Monkhorst–Pack・経路）、一般化固有値問題（正則化の規則を引数で選ぶ。既存実装の規則は、DOS・バンドは A′、jx・polB はなし）は
  [k 点と Bloch 和](kpoints-bloch.md) と [一般化固有値問題](generalized-eigen.md) のものを使います。
- **解析ごとの関数**（入出力の単位を明記し、内部は原子単位）：
  - Mulliken：\\( (\rho, S) \to M\_{I\alpha}^\sigma \\)。状態ごとの重み \\( w\_{a\mu} \\) も返せるようにする。
  - Voronoi：密度の格子と原子位置から \\( w\_I \\) で積分。\\( p \\) の合成回数（既定 6、TN9 の 3 を選択肢）と近接の範囲（近接リスト／全原子）を引数にする。
  - ESP：標本点の殻、拘束（電荷のみ／電荷と双極子）、全電荷 \\( Q\_{\mathrm{tot}} \\)（既定 0、`scf.system.charge` を選択肢）を引数にする。
  - DOS：固有値と重みの配列から、四面体法と Gauss 型の DOS・PDOS。エネルギー格子を引数で与える。
  - jx：\\( (H^\uparrow, H^\downarrow, S) \\) と対のリストから \\( J \\)。規格化（6.2 節、Cluster の \\( \tfrac12 \\) は互換オプション）、周期和／個別、極の数を引数にする。
  - Berry 位相：\\( (H, S, \mathbf D) \\) と弦の指定から \\( \varphi \\)。リンクの評価（1 次展開＋対称化／より高次）を選べるようにする。Chern 数と \\( \mathbb Z\_2 \\) は同じリンクの関数を使う。
  - 光学応答：\\( (H, S, \nabla) \\) から \\( \sigma\_{ab}(\omega) \\)、\\( \varepsilon\_{ab}(\omega) \\)。
- **既存実装の癖の扱い**：照合のための癖（Cluster の jx の \\( \tfrac12 \\)、Voronoi の 6 回、ESP の既定の倍率 \\( \lambda\_2 = 5.0 \\)、Z2FH の \\( \mathbf b \\)、calB の面積、DOS の下向きスピンの負号）は
  互換オプションで再現し、既定は物理的に正しい（または公開資料どおりの）選択にします。どちらを使ったかを出力に記録します。
- **出力形式**：`.Band`、`.Dos.val`、`.Dos.vec`、cube の形式は互換のために再現します。新しい形式（例：HDF5 や JSON）を加える場合も、単位と規約をファイル内に書きます。

## エッジケース〔規約〕〔設計〕

- **金属での Berry 位相**：占有状態の数が k 点で変わると分極は定義できません。既存実装は警告して計算を続けますが、新しい実装は既定で誤りとします。
- **縮退した固有値**：jx の周期和・Cluster の公式は \\( \varepsilon^\downarrow\_\nu = \varepsilon^\uparrow\_\mu \\) で 0/0 になります（交換分裂のない原子が縮退した準位を持つ場合など）。
  \\( \lvert\Delta\varepsilon\rvert \\) が小さい対は極限値 \\( -f'(\varepsilon) \\) で置き換えます〔設計〕。
- **非スピン偏極・非共線の入力を jx に与えた場合**：既存実装は検査せずに存在しない成分を読みます。新しい実装は誤りとします。
- **\\( T = 0 \\)**：jx と占有数の式は \\( T > 0 \\) を前提とします。
- **DOS の四面体法**：頂点のエネルギーが等しい（\\( e\_{ij} = 0 \\)）四面体と、\\( E \\) が頂点のエネルギーにちょうど一致する点では式が定義されません。
  新しい実装は縮退した四面体を、極限をとった式で扱います。
- **ESP**：van der Waals 半径の登録がない元素、殻の中に標本点がない場合、周期系（周期像の和の打ち切りの影響が大きい）。
- **Voronoi**：近接リストが自分自身だけの孤立原子（\\( w\_I = 1 \\)）。
- **`.scfout` の版とバイト順**：版 3 以外、バイト順が逆のファイル（読み手は反転して読める）、入力ファイルの節がないファイル。
- **バンドの経路**：点の数 1 の経路。

## テスト〔設計〕

`tools/check_analysis.py` は、小さな模型で次を確認しています（すべて合格）。

| 種類 | 内容 | 許容誤差 |
|---|---|---|
| 不変量 | 非直交の 1 次元鎖で、Mulliken の和則 \\( \sum\_I M\_I = N\_{\mathrm e} \\) | \\( 10^{-12} \\) |
| 不変量 | 状態ごとの Mulliken の重みの和が 1、占有数で重み付けした和が \\( \rho S \\) の総和型占有数に一致 | \\( 10^{-12} \\) |
| 不変量 | ファジーセルの重みの和が 1（合成 3 回・6 回） | \\( 10^{-14} \\) |
| 解析解 | 単純立方の s バンド（12³ の格子）：四面体法の DOS の積分が 1、\\( E\_F \\) までの積分が細かい格子での数え上げと一致 | \\( 2\times10^{-3} \\)、\\( 5\times10^{-3} \\) |
| 不変量 | Gauss 型の DOS の積分が 1 | \\( 10^{-6} \\) |
| 解析解 | Rice–Mele 鎖：格子ベクトルだけの位相と \\( e^{-i\mathbf b\cdot\mathbf r} \\) のリンクによる Berry 位相が、Resta の超格子の公式と一致（\\( \bmod 1 \\)） | \\( 10^{-10} \\) |
| 不変量 | Berry 位相のゲージ不変性（固有ベクトルに任意の位相） | \\( 10^{-12} \\) |
| 解析解 | Thouless ポンプ：1 周期の分極の変化が 1 | \\( 10^{-6} \\) |
| 解析解 | QWZ 模型の FHS Chern 数が \\( \pm1 \\)、0 | \\( 10^{-9} \\) |
| 解析解 | 3 サイトの非直交クラスター：\\( -\partial^2\Omega/\partial\theta\_i\partial\theta\_j = 2J \\)（磁気力の定理） | \\( 10^{-6} \\) |
| 不変量 | 連分数の極が Fermi 関数を再現 | \\( 10^{-12} \\) |
| 不変量 | 周期鎖で \\( \sum\_{\mathbf L}J\_{i\mathbf 0,j\mathbf L} \\)（極の和）が周期和の公式に一致 | \\( 10^{-10} \\) |

既存実装との照合（fixture、要作成）：

| 対象 | 内容 | 許容誤差の目安 |
|---|---|---|
| `.scfout` | 配布物の入力例（メタン Cluster、Si Band、bcc Fe 共線、GaAs 非共線 SOC）の `.scfout` を読み、書き直したファイルがバイト単位で一致（項目 22 の NUL 以降と、4 要素ベクトルの先頭の要素を除く） | 厳密 |
| `.scfout` | 読んだ \\( H, S, \rho \\) から \\( g\_s\sum\rho S = N\_{\mathrm e} \\)。\\( \mathbf D \\) の関係 \\( \mathbf D\_{ab}(\mathbf L) + (\mathbf R\_I - \mathbf R\_J - \mathbf L)S\_{ab}(\mathbf L) = \mathbf D\_{ba}(-\mathbf L) \\)（添字は転置） | \\( 10^{-6} \\)（格子積分の誤差） |
| Mulliken | `.out` の原子ごと・軌道ごとの値 | \\( 10^{-6} \\) |
| Voronoi | `Voronoi.charge on` の `.out`（メタン、Si） | \\( 10^{-6} \\) |
| ESP | マニュアルの例（メタン、`-c 0 -s 1.4 2.0`）：C \\( \approx -1.0696 \\)、H \\( \approx 0.2674 \\) | \\( 10^{-4} \\) |
| DOS | Si（`Dos.Kgrid 12 12 12`）の DosMain の四面体法・Gauss 型の出力 | \\( 10^{-5} \\)（相対） |
| バンド | `.Band` の固有値 | \\( 10^{-8} \\) Hartree |
| jx | bcc Fe（`Num.Kgrid` と `Num.Poles` を変えて）の \\( J \\) の周期和と個々の値 | \\( 10^{-3} \\) meV |
| polB | NaCl（マニュアル 44 節）の双極子（Debye）、Na を ±0.05 Å 動かしたときの差から \\( Z^\*\approx 1.05 \\) | \\( 10^{-6} \\) Debye |
| CDDF | Si（マニュアル 62 節）の \\( \varepsilon\_{xx}(\omega) \\) | \\( 10^{-4} \\)（相対） |

## 〔未確定〕

- jx の Cluster のモードの \\( \tfrac12 \\) が意図したものか（6.3 節）。公開資料（Han–Ozaki–Yu 2004、マニュアル Eq. 42）の係数との対応も含めて作者に確認する。
- 光学伝導度の非対角成分の符号と、帯内項の扱い（共線では入らない）、方法 1 の行列要素の詳細（非局所項の扱い）。技術ノート 2 編の式番号との照合が残っている。
- Z2FH のリンクの \\( \mathbf b \\)（8.2 節）が結果に影響する条件。
- SCF 中の Berry 位相の経路（7.3 節）とキーワードの対応。
- `DosGauss.fileout` の `.Dos.vec` の内容の詳細と、DC・NEGF の `mode 5`・`mode 6` の形式。
- 非共線の DOS の局所スピン軸への射影の式の詳細（3.2 節）。
- Band の k 経路の自動生成（マニュアル 22.1 節）、Fermi 面の出力（`FermiSurfer.fileout`、65 節）、BoltzTraP への出力（63 節）、
  バンドのアンフォールディング（55 節）、COHP、自然軌道解析（36 節）、CWF・MLWF（47 節）、XPS（60 節）、イオン化ポテンシャル（61 節）は未調査。
- `.scfout` の項目 12 と 17 が、構造最適化・MD の途中で SCF が収束しなかった場合にどの状態を書くか。
- kSpin の各モードの出力ファイルの列の詳細（9 節は概要のみ）。
