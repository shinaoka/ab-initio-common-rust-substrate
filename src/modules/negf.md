# NEGF による電子輸送

> 状態：初版（レビュー待ち）。式は技術ノート 13（Eq. 6–56）、技術ノート 4（Eq. 1–18）、技術ノート 20（Eq. 7–72）、
> T. Ozaki, PRB **75**, 035123 (2007)、T. Ozaki, K. Nishio, H. Kino, PRB **81**, 035116 (2010)、マニュアル 46 節を
> このノートの記法で書き直し、公開資料に記載のない規約・既定値は既存実装の挙動として〔規約〕で補いました。
> 極展開・表面 Green 関数・自己エネルギー・透過率・非平衡密度行列の数値的な性質は `tools/check_negf.py` で確認しています。

## この章の範囲

左右の半無限電極にはさまれた中心領域を、非平衡 Green 関数（NEGF）法で自己無撞着に解く計算（`scf.EigenvalueSolver NEGF`）と、
その結果から透過率・電流・コンダクタンス・固有チャネル・電流密度を求める後処理を扱います。

- 実空間行列と Bloch 和の規約は [H(R)・S(R) 疎ブロック形式](../design/sparse-block-format.md) と [k 点と Bloch 和](kpoints-bloch.md) に従います。
- 格子上の密度・Hartree ポテンシャル・XC は [実空間グリッド・Poisson・交換相関](grid-poisson-xc.md) を前提とし、
  この章では NEGF に固有の部分（電極からの密度の補い、境界条件付きの Poisson 方程式）だけを書きます。
- 電荷ミキシングと収束判定の一般論は [密度行列と電荷ミキシング](density-mixing.md) を参照します。
- 電極計算（Step 1）は通常の `Band` 計算そのものです。この章ではその出力の受け渡しだけを書きます。

## 記号〔設計〕

| 記号 | 意味 |
|---|---|
| \\( \mathrm L, \mathrm R \\) | 左電極、右電極。\\( s \in \\{\mathrm L, \mathrm R\\} \\) で一般の電極を表す |
| \\( \mathrm L\_0, \mathrm C\_0, \mathrm R\_0 \\) | 左電極の 1 胞分、中心部、右電極の 1 胞分。拡張中心領域 \\( \mathrm C = \mathrm L\_0\\,\vert\\,\mathrm C\_0\\,\vert\\,\mathrm R\_0 \\) |
| \\( \mathrm L\_j, \mathrm R\_j \\)（\\( j \ge 1 \\)） | C の外側に続く電極の胞。\\( \mathrm L\_1, \mathrm R\_1 \\) が C に接する |
| \\( \mathbf a\_1 \\) | 輸送方向（a 軸）。\\( \mathbf a\_2, \mathbf a\_3 \\) は面内（bc 面）の格子ベクトルで、電極と C で共通 |
| \\( \mathbf a^{(s)}\_1 \\) | 電極 \\( s \\) の単位胞の a 軸ベクトル |
| \\( \mathbf q = (q\_2, q\_3) \\) | 面内の k 点（\\( \tilde{\mathbf b}, \tilde{\mathbf c} \\) に対する分数座標） |
| \\( N\_{\mathrm C}, N\_s \\) | C の軌道数、電極 \\( s \\) の 1 胞の軌道数（非共線ではそれぞれ 2 倍） |
| \\( H\_{\mathrm C}, S\_{\mathrm C} \\) | C の \\( N\_{\mathrm C}\times N\_{\mathrm C} \\) 行列（k 点ごと） |
| \\( H^{(s)}\_{00}, S^{(s)}\_{00} \\) | 電極 \\( s \\) の 1 胞内のブロック |
| \\( H^{(s)}\_{\mathrm{out}}, S^{(s)}\_{\mathrm{out}} \\) | 電極 \\( s \\) の胞から、C から遠ざかる側の隣の胞への結合ブロック |
| \\( H\_{\mathrm Cs}, S\_{\mathrm Cs} \\) | C と \\( s\_1 \\)（\\( \mathrm L\_1 \\) または \\( \mathrm R\_1 \\)）の結合ブロック（\\( N\_{\mathrm C}\times N\_s \\)） |
| \\( z \\) | 複素エネルギー。実軸上では \\( z = E + i\eta \\)（\\( \eta > 0 \\)） |
| \\( g\_s(z) \\) | 電極 \\( s \\) の表面 Green 関数（\\( s\_1 \\) の胞のブロック） |
| \\( \Sigma\_s(z), \Gamma\_s(z) \\) | 自己エネルギー、広がり行列 |
| \\( G(z) \\) | C の遅延 Green 関数 |
| \\( \mu\_s \\) | 電極 \\( s \\) の化学ポテンシャル（バイアス込み） |
| \\( f(x) = 1/(1+e^{x}) \\)、\\( \beta = 1/(k\_{\mathrm B}T) \\) | Fermi 分布と逆温度（\\( T \\) は `scf.ElectronicTemperature`） |
| \\( V\_{\mathrm b} \\) | ソース・ドレインのバイアス（`NEGF.bias.voltage`、エネルギー単位） |

行列 \\( X \\) の Hermite 共役を \\( X^\dagger \\)、行列の「虚部」を \\( \mathrm{Im}\_{\mathrm m} X = (X - X^\dagger)/(2i) \\) と書きます。
エネルギーは内部では Hartree、入力の `NEGF.*` のエネルギーは eV です（\\( 1\\,\mathrm{Hartree} = 27.2113845\\,\mathrm{eV} \\) を使います〔規約〕）。

## 1. 系の構成

### 1.1 三領域〔公開〕

技術ノート 13（2.1 節）とマニュアル 46.1 節によれば、扱う系は a 軸方向に \\( \cdots\vert\mathrm L\_2\vert\mathrm L\_1\vert\mathrm C\vert\mathrm R\_1\vert\mathrm R\_2\vert\cdots \\) と並ぶ
一次元の鎖で、bc 面内には周期性があります。面内の Bloch 波数 \\( \mathbf q \\) ごとに独立な一次元問題になります。
局在基底が有限の半径で切れているので、次の 2 条件を満たすように胞を選べば、行列は胞についてブロック三重対角になります。

1. \\( \mathrm C\_0 \\) の基底は \\( \mathrm L\_0, \mathrm R\_0 \\) の基底とは重なってよいが、\\( \mathrm L\_1, \mathrm R\_1 \\) の基底とは重ならない。
2. 電極の各胞の基底は、隣の胞より遠い胞の基底と重ならない。

\\( \mathrm L\_0, \mathrm R\_0 \\) を C に含めるのは、電極と中心部の界面付近の電子状態の緩和を自己無撞着に取り込むためです。

### 1.2 拡張中心領域の作り方〔公開〕〔規約〕

〔公開〕マニュアル 46.3 節 A によれば、\\( \mathrm C\_0 \\) は `Atoms.Number` と `Atoms.SpeciesAndCoordinates`、
\\( \mathrm L\_0 \\) は `LeftLeadAtoms.Number` と `LeftLeadAtoms.SpeciesAndCoordinates`、
\\( \mathrm R\_0 \\) は `RightLeadAtoms.Number` と `RightLeadAtoms.SpeciesAndCoordinates` で与えます。
座標の単位は `Ang` または `AU` に限られ、`FRAC` は使えません。C の単位胞は入力しません。

〔規約〕既存実装は次のように C を組み立てます。

- **原子の並び**：C の原子は \\( \mathrm L\_0 \\)、\\( \mathrm C\_0 \\)、\\( \mathrm R\_0 \\) の順に通し番号を振ります。
- **電極との一致**：\\( \mathrm L\_0 \\) の \\( i \\) 番目の原子は、電極計算（Step 1）の単位胞の \\( i \\) 番目の原子を、
  1 番目の原子どうしが重なるように平行移動したものと、各座標 \\( 10^{-12} \\) bohr 以内で一致しなければなりません
  （\\( \mathrm R\_0 \\) も同様）。一致しなければ停止します。つまり \\( \mathrm L\_0 \\) の原子数と並びは電極の単位胞と同じです。
- **単位胞**：C の \\( \mathbf a\_2, \mathbf a\_3 \\) は左電極のものをそのまま使います（左右の電極で b・c 軸が \\( 10^{-10} \\) 以上違えば停止）。
  C の a 軸は、左電極の格子の原点を \\( \mathrm L\_0 \\) に合わせて平行移動した点 \\( \mathbf O\_{\mathrm L} \\) から、
  同様に右電極の格子の原点を \\( \mathrm R\_0 \\) に合わせた点 \\( \mathbf O\_{\mathrm R} \\) に \\( \mathbf a^{(\mathrm R)}\_1 \\) を足した点までのベクトル
  \\( \mathbf a\_1 = \mathbf O\_{\mathrm R} + \mathbf a^{(\mathrm R)}\_1 - \mathbf O\_{\mathrm L} \\) です。C の格子の原点は \\( \mathbf O\_{\mathrm L} \\) に固定されます。
  したがって C の格子の a 方向の両端の面は、左電極の胞の始まりと右電極の胞の終わりに一致します。
- **格子点数**：C の \\( N\_2, N\_3 \\) は電極の値に強制されます（左右の電極で異なれば停止）。
  \\( N\_1 \\) は通常の規則（[実空間グリッド](grid-poisson-xc.md) 1.2 節）で決まり、電極の a 方向の格子間隔とそろえる処理はありません。
- **領域の印**：標準出力に原子ごとの印を出します。`1` は \\( \mathrm C\_0 \\)、`2`／`3` は \\( \mathrm L\_0 \\)／\\( \mathrm R\_0 \\)、
  `12`／`13` は \\( \mathrm L\_0 \\)／\\( \mathrm R\_0 \\) の原子のうち \\( \mathrm L\_1 \\)／\\( \mathrm R\_1 \\) の原子と基底が重なる
  （原子間距離が両原子のカットオフ半径の和より小さい）もの。
- **重なりの検査**：\\( \mathrm C\_0 \\) の原子と \\( \mathrm L\_1, \mathrm R\_1 \\) の原子の距離がカットオフ半径の和より小さければ停止します（1.1 節の条件 1）。
  電極内の条件 2 も同様に検査します。

### 1.3 電極計算（Step 1）の条件〔公開〕〔規約〕

〔公開〕マニュアル 46.2 節：Step 1 は通常のバンド計算に `NEGF.output_hks on` と `NEGF.filename.hks` を加えたものです。
輸送方向は a 軸で、格子ベクトルの順を入れ替えて指定します（その場合 `scf.Kgrid` も入れ替える）。
半無限性を Step 2 と整合させるため、Step 1 では a 方向に十分多くの k 点を使い、bc 面の k 点数は Step 1 と Step 2 でそろえます。

〔規約〕`NEGF.output_hks on` のとき、電極の a 軸は \\( x, y, z \\) 軸のいずれかに平行で正の向きでなければなりません
（\\( (a\_x, 0, 0) \\)、\\( (0, a\_y, 0) \\)、\\( (0, 0, a\_z) \\) で成分が正）。そうでなければ停止します。
後述のゲート電圧、Poisson 方程式の線形補正、境界での \\( H \\) の混合は、この軸の Cartesian 座標を輸送方向の座標 \\( x \\) として使います。

## 2. 電極データの受け渡し

### 2.1 hks ファイルの内容〔公開〕〔規約〕

〔公開〕マニュアル 46.2 節によれば、`NEGF.filename.hks` のファイルには Hamiltonian の行列要素、電子分布、差の Hartree ポテンシャルが入ります。
Step 2 では `NEGF.filename.hks.l`、`NEGF.filename.hks.r` で左右の電極のファイルを指定します（同じファイルでもよい）。

〔規約〕既存実装が書き出す内容は、SCF 収束後の電極の次の量です（バイナリ。バイト配置は再実装で再現する必要はなく、内容だけを定めます）。

| 内容 | 用途 |
|---|---|
| 版、スピンの種類、原子数、種数、格子点数 \\( (N^{(s)}\_1, N\_2, N\_3) \\) | 整合性の検査 |
| 単位胞ベクトル、逆格子ベクトル、格子の原点、原子座標 | C の組み立て（1.2 節） |
| 化学ポテンシャル \\( \mu^{(0)}\_s \\)（Step 1 の Fermi 準位） | 電極の化学ポテンシャル |
| 種ごとの軌道数、原子ごとの隣接原子と周期像の番号 | 疎ブロック形式の構造 |
| 重なり \\( S(\mathbf L) \\)、\\( H(\mathbf L) \\)（スピンごと）、非共線では \\( H \\) の虚部（スピン軌道など） | 表面 Green 関数、自己エネルギー、境界の混合 |
| 実空間の密度行列 \\( \rho(\mathbf L) \\)（スピンごと） | 境界付近の密度の補い（5.5 節）、Mulliken 電荷の補い |
| 格子上の差の密度と差の Hartree ポテンシャル \\( \delta V\_{\mathrm H} \\) | Poisson 方程式の境界条件（7 章） |

Step 1 と Step 2 でスピンの種類（`scf.SpinPolarization`）が異なると停止します。

### 2.2 電極の k 空間ブロック〔公開〕〔規約〕

〔公開〕技術ノート 13 の Eq. 3–4 のとおり、面内の Bloch 和で k 点ごとの行列を作ります。

〔規約〕電極 \\( s \\) の実空間ブロック \\( O^{(s)}(n\_1, n\_2, n\_3) \\)（\\( \mathbf n \\) は電極の単位胞での格子ベクトルの整数座標）から

$$
O^{(s)}\_{00}(\mathbf q) = \sum\_{n\_2, n\_3} e^{2\pi i (q\_2 n\_2 + q\_3 n\_3)}\\, O^{(s)}(0, n\_2, n\_3),\qquad
O^{(s)}\_{\mathrm{out}}(\mathbf q) = \sum\_{n\_2, n\_3} e^{2\pi i (q\_2 n\_2 + q\_3 n\_3)}\\, O^{(s)}(d\_s, n\_2, n\_3)
$$

とします（\\( O = H, S \\)）。\\( d\_{\mathrm L} = -1 \\)、\\( d\_{\mathrm R} = +1 \\) で、どちらも「C から遠ざかる向き」の隣の胞への結合です。
位相の規約は [疎ブロック形式](../design/sparse-block-format.md) と同じく格子ベクトルだけに付けます。逆向きの結合は
\\( O^{(s)}\_{\mathrm{in}} = (O^{(s)}\_{\mathrm{out}})^\dagger \\) です。

C と \\( s\_1 \\) の結合 \\( H\_{\mathrm Cs}, S\_{\mathrm Cs} \\) も電極のデータから作ります。\\( \mathrm L\_0 \\) は電極の胞と同一なので、
\\( H\_{\mathrm C\mathrm L} \\) の \\( \mathrm L\_0 \\) の行には \\( H^{(\mathrm L)}\_{\mathrm{out}} \\) をそのまま置き、\\( \mathrm C\_0, \mathrm R\_0 \\) の行は 0 とします
（右も同様）。C の側で計算した行列要素は、この結合には使いません。

### 2.3 C の行列と境界での混合〔規約〕

C の \\( H\_{\mathrm C}(\mathbf q), S\_{\mathrm C}(\mathbf q) \\) は、C の SCF で得た実空間行列のうち **a 方向の周期像を含まない**項
（\\( n\_1 = 0 \\)）だけを Bloch 和したものです。a 方向の隣の像との結合は、自己エネルギーで置き換えられます。

さらに、印 `12` の原子どうし（\\( \mathrm L\_0 \\) のうち \\( \mathrm L\_1 \\) と重なる原子）のブロックは

$$
S\_{\mathrm C} \leftarrow S^{(\mathrm L)}\_{00},\qquad
H\_{\mathrm C} \leftarrow c\\, H\_{\mathrm C} + (1 - c)\\, H^{(\mathrm L)}\_{00}
$$

と電極の値で置き換え・混合します（`13` と右電極も同様）。重み \\( c \\) は `NEGF.Mixing.H.Boundary` が 0 以上ならその定数、
既定（\\( -1 \\)）では原子対の輸送方向座標の中点 \\( \bar x \\) と C の端 \\( x\_{\mathrm e} \\)（左は格子の原点、右は原点 \\( + N\_1 h\_1 \\)）から

$$
c(\bar x) = \tfrac12\left[\tanh\\!\left(\tfrac12\left(\lvert\bar x - x\_{\mathrm e}\rvert - 2\right)\right) + 1\right]\qquad(\bar x\\ \text{は bohr})
$$

です。端から 2 bohr の位置で \\( c = 1/2 \\)、端では \\( c \approx 0.12 \\)、端から 8 bohr で \\( c \approx 0.998 \\) です。
電極の \\( H^{(s)}\_{00} \\) はバイアスを加えた後の値（6.1 節）です。

## 3. 表面 Green 関数〔公開〕〔規約〕

〔公開〕技術ノート 13（2.2 節、Eq. 9–21）は、表面 Green 関数 \\( g\_s(z) \\) を繰り込み（decimation）法で求めます。
電極 \\( s \\) について、最初の値を

$$
\varepsilon\_{\mathrm s}^{(0)} = \varepsilon^{(0)} = zS^{(s)}\_{00} - H^{(s)}\_{00},\qquad
\alpha^{(0)} = -\left(zS^{(s)}\_{\mathrm{out}} - H^{(s)}\_{\mathrm{out}}\right),\qquad
\beta^{(0)} = -\left(zS^{(s)}\_{\mathrm{in}} - H^{(s)}\_{\mathrm{in}}\right)
$$

とし、\\( i = 0, 1, \dots \\) について

$$
\begin{aligned}
\varepsilon\_{\mathrm s}^{(i+1)} &= \varepsilon\_{\mathrm s}^{(i)} - \alpha^{(i)} (\varepsilon^{(i)})^{-1} \beta^{(i)}, &
\varepsilon^{(i+1)} &= \varepsilon^{(i)} - \alpha^{(i)} (\varepsilon^{(i)})^{-1} \beta^{(i)} - \beta^{(i)} (\varepsilon^{(i)})^{-1} \alpha^{(i)},\\\\
\alpha^{(i+1)} &= \alpha^{(i)} (\varepsilon^{(i)})^{-1} \alpha^{(i)}, &
\beta^{(i+1)} &= \beta^{(i)} (\varepsilon^{(i)})^{-1} \beta^{(i)}
\end{aligned}
$$

を繰り返すと \\( g\_s = \lim\_{i\to\infty} (\varepsilon\_{\mathrm s}^{(i)})^{-1} \\) です。\\( i \\) 回目の繰り込みは \\( 2^i \\) 個の胞を消去したことに当たり、
\\( \lVert\alpha^{(i)}\rVert, \lVert\beta^{(i)}\rVert \\) は \\( \mathrm{Im}\\,z > 0 \\) なら二重指数的に 0 に近づきます。
得られた \\( g\_s \\) は Dyson 方程式

$$
g\_s = \left[zS^{(s)}\_{00} - H^{(s)}\_{00} - \left(zS^{(s)}\_{\mathrm{out}} - H^{(s)}\_{\mathrm{out}}\right) g\_s \left(zS^{(s)}\_{\mathrm{in}} - H^{(s)}\_{\mathrm{in}}\right)\right]^{-1}
$$

を満たします（`tools/check_negf.py` で残差 \\( < 10^{-13} \\) を確認）。

〔規約〕

- 収束判定は、連続する 2 回の \\( (\varepsilon\_{\mathrm s}^{(i)})^{-1} \\) の要素ごとの差の絶対値の最大が
  `NEGF.Surfgreen.convergeeps`（既定 \\( 10^{-12} \\)）未満になったときです。最大反復回数は `NEGF.Surfgreen.iterationmax`（既定 600）で、
  到達しても停止せず、エラーの行を出力して最後の値を使います。
- \\( z \\) は密度行列・透過率のどの段階でも \\( \mathrm{Im}\\,z > 0 \\) か \\( \mathrm{Im}\\,z < 0 \\) のどちらかで、実軸上で直接は評価しません。
  先進側（\\( \mathrm{Im}\\,z < 0 \\)）が必要な場合は、\\( z^\* \\) で同じ繰り込みを行います（2 回目の計算をする）。
- 繰り込み法のほかに、固有値問題による方法（Sanvito らの方法）の実装も残っていますが、使われていません。

## 4. 自己エネルギーと中心の Green 関数〔公開〕

技術ノート 13 の Eq. 6–8 を書き直すと、

$$
\Sigma\_s(z) = \left(zS\_{\mathrm Cs} - H\_{\mathrm Cs}\right) g\_s(z) \left(zS\_{s\mathrm C} - H\_{s\mathrm C}\right),\qquad
S\_{s\mathrm C} = S\_{\mathrm Cs}^\dagger,\\ H\_{s\mathrm C} = H\_{\mathrm Cs}^\dagger
$$

$$
G(z) = \left[zS\_{\mathrm C} - H\_{\mathrm C} - \Sigma\_{\mathrm L}(z) - \Sigma\_{\mathrm R}(z)\right]^{-1}
$$

です。非直交基底なので、自己エネルギーは結合の重なり \\( S\_{\mathrm Cs} \\) を通じて \\( z \\) に陽に依存します。
\\( G(z) \\) は \\( z \\) について上半面・下半面で解析的で、\\( G(z^\*) = G(z)^\dagger \\) です。広がり行列は

$$
\Gamma\_s(E) = i\left[\Sigma\_s(E + i\eta) - \Sigma\_s(E + i\eta)^\dagger\right]
$$

です（技術ノート 13 Eq. 45、技術ノート 20 Eq. 7）。

〔規約〕\\( \Sigma\_s \\)、\\( G \\) は密行列として作り、\\( G \\) は LU 分解による逆行列で求めます。共線スピン偏極ではスピンごと、
非共線では \\( 2N \\) 次元（12 章）で計算します。

## 5. 密度行列

### 5.1 平衡部分：Fermi 関数の極展開〔公開〕

技術ノート 4（Eq. 1–10）と PRB **75**, 035123 は、Fermi 関数の連分数表示を有限の深さ \\( M = 2N\_{\mathrm p} \\) で打ち切り、部分分数に分解します。
このノートの記法で導出を書きます。打ち切った連分数は

$$
f(x) \simeq \frac12 - \frac{x}{4}\left[(D + ixB)^{-1}\right]\_{11},\qquad
D = \mathrm{diag}(1, 3, 5, \dots, 2M - 1),\quad B\_{p,p\pm 1} = \tfrac12
$$

と書けます（\\( B \\) は対角が 0 の三重対角）。一般化固有値問題 \\( Bv\_p = \lambda\_p D v\_p \\)、\\( v\_p^{\mathsf T}Dv\_p = 1 \\) を解くと
\\( (D + ixB)^{-1} = \sum\_p v\_p v\_p^{\mathsf T}/(1 + ix\lambda\_p) \\) で、固有値は \\( \pm\lambda \\) の対で現れ、対の \\( v\_{1p}^2 \\) は等しくなります。
対をまとめると

$$
f(x) \simeq \frac12 + \sum\_{p=1}^{N\_{\mathrm p}} R\_p\left[\frac{1}{x - iz\_p} + \frac{1}{x + iz\_p}\right],\qquad
z\_p = \frac{1}{\lambda\_p}\\ (\lambda\_p > 0),\quad R\_p = -\frac{v\_{1p}^2\\, z\_p^2}{4}
$$

です（技術ノート 4 Eq. 3）。\\( z\_p \\) は実数で、小さい方の約 6 割は Matsubara 極 \\( (2p - 1)\pi \\) に、\\( R\_p \\) は \\( -1 \\) に一致し、
残りの極は急速に疎になって遠方に伸びます（`tools/check_negf.py`：\\( N\_{\mathrm p} = 20 \\) で既に \\( \lvert x\rvert \le 60 \\) の誤差 \\( 10^{-15} \\)）。

\\( x = \beta(E - \mu) \\) とし、\\( \alpha\_p = \mu + iz\_p/\beta \\) とおくと、\\( G(z) \\) の極が実軸上（または実軸の近く）にある場合、

$$
\rho\_{\mathrm{eq}}(\mu) = -\frac{1}{\pi}\mathrm{Im}\_{\mathrm m}\\!\int\_{-\infty}^{\infty} dE\\, G(E + i0^+) f\\big(\beta(E - \mu)\big)
{}= \frac12\\, m^{(0)} - \frac{1}{\beta}\sum\_{p=1}^{N\_{\mathrm p}} R\_p\left[G(\alpha\_p) + G(\alpha\_p)^\dagger\right]
$$

です（スピンあたり）。\\( R\_p < 0 \\) なので各項の重み \\( -R\_p/\beta \\) は正です。
技術ノート 4 Eq. 13 はスピンの因子 2 を含めた形（\\( \mu^{(0)} + \mathrm{Im}[-4i\beta^{-1}\sum G(\alpha\_p)R\_p] \\)）で、上の式の 2 倍です。
孤立系（自己エネルギーなし）では \\( \rho\_{\mathrm{eq}} = \sum\_n c\_n f(\varepsilon\_n) c\_n^\dagger \\) と一致します（`tools/check_negf.py` で \\( 10^{-15} \\)）。

### 5.2 0 次モーメント〔公開〕〔規約〕

〔公開〕\\( m^{(0)} = -\pi^{-1}\mathrm{Im}\_{\mathrm m}\int G(E+i0^+)\\,dE \\) は Green 関数の 0 次モーメントで、
大きな \\( R \\) について \\( m^{(0)} \simeq iR\\,G(iR) \\) と評価します（技術ノート 4 Eq. 14–18、技術ノート 13 Eq. 26）。

〔規約〕既存実装は \\( R = 10^{12} \\) Hartree とし、\\( G(iR) \\) を自己エネルギー込みで（表面 Green 関数も \\( z = iR \\) で求めて）計算します。
\\( S\_{\mathrm Cs} \ne 0 \\) のとき \\( \Sigma\_s(z) \\) は \\( z \\) に比例して大きくなるので、\\( m^{(0)} \\) は \\( S\_{\mathrm C}^{-1} \\) とは一致しません
（`tools/check_negf.py` の模型で差 \\( 7\times 10^{-3} \\)）。再実装でも \\( S\_{\mathrm C}^{-1} \\) で代用してはいけません。

〔規約〕平衡部分で \\( G \\) を評価する点の数は \\( N\_{\mathrm p} + 1 \\) です（\\( N\_{\mathrm p} \\) は `NEGF.Num.Poles`、既定 150）。

### 5.3 有限バイアス：非平衡部分〔公開〕〔規約〕

〔公開〕技術ノート 13（3.1 節、Eq. 42–46）とマニュアル 46.3 節 B によれば、バイアス下の密度行列は、一方の電極の化学ポテンシャルでの
平衡部分と、実軸上の非平衡の補正の和です。非平衡部分は \\( z \\) について解析的でないので、
実軸から小さな虚部 \\( \eta\_{\mathrm{neq}} \\) だけ離した線上で、等間隔の矩形則で積分します。

〔規約〕既存実装の組み立ては次のとおりです。\\( \mu\_< = \min(\mu\_{\mathrm L}, \mu\_{\mathrm R}) \\)、\\( \mu\_> = \max(\mu\_{\mathrm L}, \mu\_{\mathrm R}) \\) とし、
\\( \mu\_> \\) を持つ電極を \\( s\_> \\) とします。

$$
\rho = \rho\_{\mathrm{eq}}(\mu\_<) + \frac{\Delta E}{2\pi}\sum\_{i=0}^{n\_{\mathrm{neq}}-1}
\left[f\\big(\beta(E\_i - \mu\_>)\big) - f\\big(\beta(E\_i - \mu\_<)\big)\right] G(z\_i)\\,\Gamma\_{s\_>}(z\_i)\\,G(z\_i)^\dagger,\qquad z\_i = E\_i + i\eta\_{\mathrm{neq}}
$$

ここで \\( \Gamma\_{s\_>}(z) = i[\Sigma\_{s\_>}(z) - \Sigma\_{s\_>}(z)^\dagger] \\)、
\\( \Delta E \\) は `NEGF.bias.neq.energy.step`（既定 0.02 eV）、\\( \eta\_{\mathrm{neq}} \\) は `NEGF.bias.neq.im.energy`（既定 0.01 eV）です。
平衡部分の極は \\( \mu\_< \\) を中心に置きます。

積分区間は、\\( \bar\mu = (\mu\_< + \mu\_>)/2 \\) について対称な \\( [E\_{\min}, E\_{\max}] \\) です。\\( E\_{\max} \\) は、
\\( \bar\mu \\) から上へ \\( \Delta E/10 \\) 刻みの格子上で、\\( \lvert f(\beta(E - \mu\_>)) - f(\beta(E - \mu\_<))\rvert \\) が
`NEGF.bias.neq.cutoff`（既定 \\( 10^{-8} \\)）未満になる最初の \\( E \\) とし、\\( E\_{\min} = 2\bar\mu - E\_{\max} \\) です。
点数は \\( n\_{\mathrm{neq}} = \lfloor (E\_{\max} - E\_{\min})/\Delta E\rfloor \\)、点は \\( E\_i = E\_{\min} + i\\,\Delta E \\)（\\( i = 0, \dots, n\_{\mathrm{neq}} - 1 \\)）です。

\\( V\_{\mathrm b} = 0 \\) では最初の点で差が 0 なので \\( n\_{\mathrm{neq}} = 0 \\) となり、平衡部分だけになります。
\\( G \\) を評価する点の総数は \\( N\_{\mathrm p} + 1 + n\_{\mathrm{neq}} \\) で、計算時間はこれに比例します（マニュアル 46.3 節）。
例：マニュアル 46.3 節の出力例（炭素鎖、\\( V\_{\mathrm b} = 0.5 \\) V）の区間は \\( E\_{\max} - \mu\_> = 0.954 \\) eV、\\( n\_{\mathrm{neq}} = 120 \\) で、例題入力の電子温度 600 K の \\( 18.4\\,k\_{\mathrm B}T \\)（\\( e^{-18.4} \approx 10^{-8} \\)）に当たります。

〔設計〕参照する電極の選び方は任意で、厳密には \\( \rho\_{\mathrm{eq}}(\mu\_>) - \int G\Gamma\_{s\_<}G^\dagger\\,\Delta f\\,dE/2\pi \\) と同じです
（\\( G(\Gamma\_{\mathrm L} + \Gamma\_{\mathrm R})G^\dagger = i(G - G^\dagger) \\) で、束縛状態がない場合）。
`tools/check_negf.py` ではこの 2 通りが \\( \eta\_{\mathrm{neq}} \\) と刻みの誤差の範囲（\\( 10^{-4} \\)）で一致することを確認しています。
\\( \eta\_{\mathrm{neq}} \\) は区間内の束縛状態の寄与を \\( \eta \\) 幅のローレンツ関数で拾うための正則化でもあり、\\( \Delta E \lesssim \eta\_{\mathrm{neq}} \\) 程度でないと積分が粗くなります。

### 5.4 k 点と実空間の密度行列〔公開〕〔規約〕

〔公開〕技術ノート 13 Eq. 22, 43 は、面内の Brillouin 帯域の積分で実空間の密度行列を作ります（同ノートは位相 \\( e^{-i\mathbf k\cdot\mathbf R} \\) で書いています）。
マニュアル 46.3 節：面内の k 点数は `NEGF.scf.Kgrid`（既定 `1 1`）で、a 方向には指定しません。

〔規約〕

- k 点は \\( q\_j = -\tfrac12 + \tfrac{2l + 1}{2N\_j} \\)（\\( j = 2, 3 \\)、\\( l = 0, \dots, N\_j - 1 \\)）で、[k 点](kpoints-bloch.md) の `Regular` と同じ式です
  （\\( q\_2 \\) に \\( +10^{-12} \\)、\\( q\_3 \\) に \\( -10^{-12} \\) の微小なずれを加える）。Monkhorst–Pack の縮約（`scf.MP.criterion` の方式）は NEGF では使えず、指定すると停止します。
- 共線計算では \\( (l\_2, l\_3) \\) と \\( (N\_2 - 1 - l\_2, N\_3 - 1 - l\_3) \\) を時間反転の対として一方だけを重み 2 で計算し、自分自身が相手の点は重み 1 です。
  非共線計算では削減しません。
- 実空間の密度行列は、k 点ごとの \\( \rho(\mathbf q) \\)（5.1–5.3 節の行列）から
  \\( \rho\_{I\alpha,J\beta}(\mathbf n) = \frac{1}{N\_2N\_3}\sum\_{\mathbf q} w\_{\mathbf q}\\,\mathrm{Re}\left[\rho\_{I\alpha,J\beta}(\mathbf q)\\,e^{-2\pi i(q\_2n\_2 + q\_3n\_3)}\right] \\)
  とします（虚部も同じ位相で別に保持）。\\( \rho(\mathbf q)\_{I\alpha,J\beta} \\) は \\( c\_{I\alpha}\overline{c\_{J\beta}} \\) の形なので、
  これは [疎ブロック形式](../design/sparse-block-format.md) の \\( e^{+i\mathbf k\cdot\mathbf L}\\,\overline{c\_{I\alpha}}c\_{J\beta} \\) と同じ規約です。
- C の原子と a 方向の周期像（\\( n\_1 = \pm 1 \\)）の間の密度行列要素は 0 にします。
- NEGF の段階では、化学ポテンシャルの探索（電子数の拘束）はしません。C の電子数は電極の化学ポテンシャルで決まる結果です。

### 5.5 格子上の密度と電極からの補い〔公開〕〔規約〕

〔公開〕技術ノート 13（2.4 節、Eq. 27–31）によれば、C の格子上の密度は、C の基底どうしの項、C と電極の基底の交差項、電極の基底どうしの項の和です。

〔規約〕既存実装は次のように作ります。

- C の密度行列（5.4 節）からの寄与は、通常の格子上の密度（[実空間グリッド](grid-poisson-xc.md) 2.1 節）と同じですが、
  軌道の a 方向の周期像（\\( n\_1 \ne 0 \\)）からの寄与を含めません。原子密度の重ね合わせ（\\( \delta n \\) の基準）も同様です。
- 両端の境界の板状の領域には、\\( \mathrm L\_1 \\)（\\( \mathrm R\_1 \\)）の原子の基底の裾の寄与を **Step 1 の電極の密度行列** で加えます：
  \\( \mathrm L\_1 \\) の軌道どうしの項と、\\( \mathrm L\_1 \\) と \\( \mathrm L\_0 \\) の軌道の交差項の 2 倍です。
  面内は \\( \mathbf a\_2, \mathbf a\_3 \\) 方向に \\( -1, 0, +1 \\) の像まで足します。
  境界の領域は、左は a 方向の格子面 0 から「\\( \mathrm L\_0 \\) の原子のカットオフ球の最も右の端を \\( -\mathbf a^{(\mathrm L)}\_1 \\) だけ動かした面 \\( + 1 \\)」まで、右も対称に定めます。
  原子密度の重ね合わせにも、\\( \mathrm L\_1, \mathrm R\_1 \\) の原子の寄与を同じ領域に加えます。
- Mulliken 電荷には、\\( \mathrm L\_0 \\)（\\( \mathrm R\_0 \\)）の原子の軌道について、\\( \mathrm L\_1 \\)（\\( \mathrm R\_1 \\)）の軌道との
  \\( \sum\_\beta S\_{\alpha\beta}\rho\_{\alpha\beta} \\) を Step 1 の電極の \\( S \\)、\\( \rho \\) から加えます。

### 5.6 旧積分経路（`NEGF.Integration OLD`）〔規約〕

既定（`CF`）は 5.1 節の極展開です。`OLD` を指定すると、\\( \mu \\) から虚部 \\( \Delta = 1.2\\,\mathrm{eV} \\) 以内の Matsubara 極、
虚部 \\( \Delta \\) の線分（\\( \mu \pm 40k\_{\mathrm B}T \\)、Gauss–Legendre 400 点）、下端 \\( -1 \\) Hartree に戻る円弧（40 点）の経路を使います。
バイアス 0 でしか使えず、指定すると有限バイアスでは停止します。技術ノート 13 付録の経路（\\( \Delta = 1.2 \\) eV、\\( \gamma = 40k\_{\mathrm B}T \\)、40 点の円弧）はこれに当たります。
再実装では対象外とします〔設計〕。

## 6. バイアスとゲート

### 6.1 バイアス電圧〔公開〕〔規約〕

〔公開〕技術ノート 13 Eq. 47–49：右電極のバイアスは \\( H\_{\mathrm R} \to H\_{\mathrm R} + V\_{\mathrm b}S\_{\mathrm R} \\)、\\( \mu\_{\mathrm R} \to \mu\_{\mathrm R} + V\_{\mathrm b} \\)、
境界の差の Hartree ポテンシャルを同じだけずらすことで与えます。電極の静電ポテンシャルの一様なずれは、局在基底では \\( V S \\) の加算と等価です。
マニュアル 46.3 節：`NEGF.bias.voltage` は eV で与え、1 eV の化学ポテンシャル差を 1 V と読みます。

〔規約〕既存実装は電圧を左右に半分ずつ振り分けます。

$$
V\_{\mathrm L} = -\tfrac12 V\_{\mathrm b},\quad V\_{\mathrm R} = +\tfrac12 V\_{\mathrm b};\qquad
H^{(s)}(\mathbf L) \to H^{(s)}(\mathbf L) + V\_s S^{(s)}(\mathbf L),\quad
\mu\_s = \mu^{(0)}\_s + V\_s,\quad
\delta V^{(s)}\_{\mathrm H}(\mathbf r) \to \delta V^{(s)}\_{\mathrm H}(\mathbf r) + V\_s
$$

（電子のポテンシャルエネルギーとしての加算。正のバイアスで右電極の準位が上がる）。共線ではすべてのスピン成分、
非共線では \\( \uparrow\uparrow, \downarrow\downarrow \\) の成分にだけ加えます（\\( V\_s\\,I\_2\otimes S \\)）。
加える対象は電極の全ブロック（\\( H\_{00} \\)、\\( H\_{\mathrm{out}} \\)、C との結合）です。

> **マニュアル 46.3 節についての注記**：同節の出力例は「左に 0、右に \\( V\_{\mathrm b} \\)」を加えた形ですが、参照した実装は \\( \pm V\_{\mathrm b}/2 \\) です。
> \\( \mu\_{\mathrm R} - \mu\_{\mathrm L} \\) と非平衡の窓の幅は同じで、全体のエネルギーの原点が \\( V\_{\mathrm b}/2 \\) ずれるだけです。
> ただし C の静電ポテンシャルの基準（7 章の境界値）がこれに従って変わるので、照合では振り分け方をそろえる必要があります。

`NEGF.bias.apply`（既定 `on`）を `off` にするとバイアスの処理全体を省きます（この経路には末尾の〔未確定〕に挙げた疑問があります）。

### 6.2 ゲート電圧〔公開〕〔規約〕

〔公開〕技術ノート 13 Eq. 50、マニュアル 46.3 節 D：

$$
V\_{\mathrm g}(x) = V\_{\mathrm g}^{(0)}\exp\\!\left[-\left(\frac{x - x\_{\mathrm c}}{d}\right)^{8}\right]
$$

を加えます。\\( V\_{\mathrm g}^{(0)} \\) は `NEGF.gate.voltage`（eV、既定 0）です。

〔規約〕\\( x \\) は格子面の輸送方向の Cartesian 座標（面の番号 \\( \times\\,\lvert\mathbf h\_1\rvert \\) + 原点）、
\\( x\_{\mathrm c} \\) は \\( \mathrm C\_0 \\) の原子の \\( x \\) 座標の平均、\\( d = \lvert a\_{1,x}\rvert - \lvert a^{(\mathrm L)}\_{1,x}\rvert - \lvert a^{(\mathrm R)}\_{1,x}\rvert \\)
（\\( \mathrm C\_0 \\) の部分の長さ）です。\\( V\_{\mathrm g} \\) は Poisson 方程式を解いた後の差の Hartree ポテンシャル（電子のポテンシャルエネルギー）に、
SCF の各反復で加えます。したがって \\( V\_{\mathrm g}^{(0)} > 0 \\) は \\( \mathrm C\_0 \\) の準位を押し上げます。

### 6.3 非 SCF のバイアス（`NEGF.NonSCF.bias`）〔規約〕

`on`（既定 `off`）のとき、C の \\( H\_{\mathrm C} \\) に、原子対の中点 \\( \bar x \\) での一次関数
\\( V(\bar x) = V\_{\mathrm L} + (V\_{\mathrm R} - V\_{\mathrm L})(\bar x - x\_0)/(x\_1 - x\_0) \\) を \\( V(\bar x)\\,S\_{\mathrm C} \\) として加えます（\\( x\_0, x\_1 \\) は C の格子の両端）。

## 7. Poisson 方程式の境界条件

### 7.1 二次元 FFT ＋ 一次元差分（`FD`、既定）〔公開〕〔規約〕

〔公開〕技術ノート 13（2.5 節、Eq. 32–37）：差の密度 \\( \delta n \\) に対する Poisson 方程式を、bc 面内で Fourier 変換し、
a 方向は差分で解きます。面内の逆格子ベクトルを \\( \mathbf G\_\parallel \\) とすると

$$
\left(\frac{d^2}{dx^2} - \lvert\mathbf G\_\parallel\rvert^2\right)\delta\tilde V\_{\mathrm H}(x, \mathbf G\_\parallel) = -4\pi\\,\delta\tilde n(x, \mathbf G\_\parallel)
$$

で、両端の値を電極の Step 1 の \\( \delta V\_{\mathrm H} \\) に固定する Dirichlet 条件です。界面での静電ポテンシャルが電極と連続になるので、
理論的にはこちらが適切です（マニュアル 46.3 節）。

〔規約〕

- 格子面 \\( k = 0, \dots, N\_1 - 1 \\)、間隔 \\( h\_1 = \lvert\mathbf h\_1\rvert \\) で
  \\( V\_{k-1} - (2 + h\_1^2\lvert\mathbf G\_\parallel\rvert^2)V\_k + V\_{k+1} = -4\pi h\_1^2\\,\delta\tilde n\_k \\)
  を、\\( \mathbf G\_\parallel \\) ごとに三重対角の連立一次方程式として解きます（\\( h\_1 \\) は a 軸の格子ベクトルの長さで、a 軸が bc 面に直交することを暗に仮定）。
- 境界値 \\( V\_{-1} \\) は左電極の \\( \delta V\_{\mathrm H} \\)（バイアス \\( V\_{\mathrm L} \\) 込み）の **最後の** a 面、\\( V\_{N\_1} \\) は右電極の **最初の** a 面を
  bc 面で Fourier 変換した値です。C の格子の両端が電極の胞の境界に一致する（1.2 節）ので、これらは C の外側に隣接する面です。
- `NEGF.Poisson_Gparazero.scaling`（既定 1.0）は \\( \mathbf G\_\parallel = 0 \\) の成分の \\( \delta\tilde n \\) に掛ける係数です。
- 解いた後、6.2 節のゲート電圧を加えます。
- 電荷ミキシングは、この表現（bc 面の Fourier 成分と a 方向の格子面）の上で行う NEGF 用の変種を使います（詳細は〔未確定〕）。

### 7.2 三次元 FFT（`FFT`）〔公開〕〔規約〕

〔公開〕マニュアル 46.3 節：三次元 FFT（周期境界）で解く方法も選べます。非極性の系では FD とほぼ同じ結果で、SCF は収束しやすい傾向があります。

〔規約〕通常の三次元 FFT の解（\\( \mathbf G = 0 \\) 成分は 0）に、a 方向の一次関数 \\( v(x) = v\_0 + (v\_1 - v\_0)(x - x\_0)/(x\_1 - x\_0) \\) を加えます。
\\( v\_0 \\)（\\( v\_1 \\)）は、左（右）電極の \\( \delta V\_{\mathrm H} \\) の最初の a 面の面平均から、FFT 解の C の最初（最後）の a 面の面平均を引いた値で、
\\( x\_0, x\_1 \\) は C の格子の両端の座標です。境界で一致するのは面平均（\\( \mathbf G\_\parallel = 0 \\) 成分）だけです。

その他、`FDG`（二次元の差分とゲート用の 2 つの境界値、`NEGF.gate.voltage` に 2 つの値を与える）という未記載の方式があります〔規約〕。再実装では対象外とします〔設計〕。

### 7.3 Band から NEGF への切り替え時の Hartree ポテンシャル〔規約〕

`NEGF.Band2NEGF.Hartree on`（既定 `off`）のとき、FD の場合に限り、NEGF に切り替えた後の最初の 5 反復で、
解いた \\( \delta\tilde V\_{\mathrm H} \\) を、両境界値の間の一次補間（面 \\( -1 \\) と \\( N\_1 \\) を結ぶ直線）と混ぜます。
解の重みは 1–5 回目でそれぞれ 0.05、0.15、0.25、0.50、0.75 です。

## 8. SCF の流れ〔公開〕〔規約〕

〔公開〕マニュアル 46.1・46.3 節によれば、(i) 電極の SCF（Step 1）、(ii) C の平衡 Green 関数、有限バイアスでは (iii) 非平衡 Green 関数、
(iv) 密度の更新と混合、を収束まで繰り返し、その後に透過率などを計算します（Step 3）。
最初の数回は a 方向にも周期性を仮定した通常の対角化（Band）で回すと、最終的な収束が速くなります（`NEGF.scf.Iter.Band`）。
NEGF の収束判定は `scf.criterion` を残差ノルム `NormRD` に適用します（他の解法では `dUele`）。

〔規約〕

- SCF の反復番号が `NEGF.scf.Iter.Band` 以下の間は Band 計算（a 方向は Γ 点だけ、面内は `NEGF.scf.Kgrid`）で密度行列を作り、
  それ以降は NEGF（5 章）で作ります。既定値は **10** です（マニュアルの 6 と異なる）。密度行列ファイルから再開した場合は最初から NEGF です。
- Band の段階の化学ポテンシャルは、`NEGF.ChemP.Band fix`（既定）なら \\( (\mu\_{\mathrm L} + \mu\_{\mathrm R})/2 \\) に固定、`free` なら電子数から決めます。
- Poisson 方程式は Band の段階から 7 章の方式で解きます。
- 収束判定は残差ノルムの平方根 \\( \sqrt{\lvert\mathrm{NormRD}\rvert} \\) が `scf.criterion` 未満のときです。
- NEGF の段階ではバンドエネルギーの項（\\( E\_{\mathrm{band}} \\) に当たる量）を 0 とし、エネルギー重み付き密度行列も計算しません（9.5 節）。

## 9. 透過率・電流・コンダクタンス（Step 3）

### 9.1 透過率〔公開〕〔規約〕

〔公開〕技術ノート 13 Eq. 51–52（Caroli の式）：

$$
T\_\sigma(E; \mathbf q) = \mathrm{Tr}\left[\Gamma\_{\mathrm L}(E)\\,G(E + i\eta)\\,\Gamma\_{\mathrm R}(E)\\,G(E - i\eta)\right],\qquad
T\_\sigma(E) = \frac{1}{N\_{\mathbf q}}\sum\_{\mathbf q} T\_\sigma(E; \mathbf q)
$$

k 点は `NEGF.tran.Kgrid`（既定は `NEGF.scf.Kgrid` の値）で、SCF より多くすることが推奨されています（マニュアル 46.4.1 節）。

〔規約〕

- エネルギーの格子は \\( E\_i = \mu\_{\mathrm L} + E\_{\mathrm{lo}} + (E\_{\mathrm{hi}} - E\_{\mathrm{lo}})\\,i/(n\_{\mathrm{div}} - 1) \\)（\\( i = 0, \dots, n\_{\mathrm{div}} - 1 \\)）で、
  `NEGF.tran.energyrange` の 3 つの値が \\( E\_{\mathrm{lo}}, E\_{\mathrm{hi}}, \eta \\)（既定 \\( -10, 10, 10^{-3} \\) eV）、\\( n\_{\mathrm{div}} \\) は `NEGF.tran.energydiv`（既定 200）です。
  エネルギーの原点はバイアス込みの \\( \mu\_{\mathrm L} \\) です。
- \\( G(E - i\eta) \\) と \\( \Sigma\_s(E - i\eta) \\) は \\( z^\* \\) で計算し直し、\\( \Gamma\_s = i[\Sigma\_s(E + i\eta) - \Sigma\_s(E - i\eta)] \\) とします（数学的には \\( \Sigma^\dagger \\) と同じ）。トレースの実部を出力します。
- k 点の格子と時間反転による削減は 5.4 節と同じで、削減した点の値を相手の点にも写します。
- 有限の \\( \eta \\) は C の中で吸収として働き、完全な鎖でも \\( 1 - T \\) は \\( \eta \\) に比例して増えます
  （`tools/check_negf.py` の 5 胞の鎖で \\( \eta = 10^{-5} \\) Hartree のとき \\( 1 - T = 4\times 10^{-4} \\)）。既定の \\( 10^{-3} \\) eV（\\( 3.7\times 10^{-5} \\) Hartree）での照合では、この系統誤差を見込みます。

### 9.2 電流〔公開〕〔規約〕

〔公開〕技術ノート 13 Eq. 53（Landauer–Büttiker）：

$$
I\_\sigma = \frac{e}{h}\int dE\\, T\_\sigma(E)\left[f\\big(\beta(E - \mu\_{\mathrm R})\big) - f\\big(\beta(E - \mu\_{\mathrm L})\big)\right]
$$

（スピンあたり。符号は \\( \mu\_{\mathrm R} > \mu\_{\mathrm L} \\) で正）。

〔規約〕

- 原子単位（\\( e = \hbar = 1 \\)、\\( e/h = 1/2\pi \\)）で \\( I\_\sigma = -\frac{\Delta E}{2\pi}\sum\_i [f\_{\mathrm L}(E\_i) - f\_{\mathrm R}(E\_i)]\\,T\_\sigma(E\_i + i\eta\_{\mathrm I}) \\) を計算し、
  \\( eE\_{\mathrm h}/\hbar = 6.6236178\times 10^{-3} \\) A を掛けてアンペアにします（CODATA 2018 では \\( 6.6236182\times 10^{-3} \\) A。相対差 \\( 7\times 10^{-8} \\)）。
- 積分の区間と点は 5.3 節と同じ規則で、刻み \\( \Delta E \\) は `NEGF.current.energy.step`（既定 0.01 eV）、打ち切りは `NEGF.current.cutoff`（既定 \\( 10^{-8} \\)）、
  虚部 \\( \eta\_{\mathrm I} \\) は `NEGF.current.im.energy`（既定は `NEGF.tran.energyrange` の 3 番目の値）です。点数が 10 未満なら刻みが粗いという警告を出します。
- 温度は `scf.ElectronicTemperature`（既定 300 K）。
- 電流は k 点ごとに計算し、`System.Name.current` に k 点ごとの値と平均を出します。スピン偏極なしでは両スピンの欄に同じ値（スピンあたりの値）を出すので、全電流はその 2 倍です。

### 9.3 コンダクタンス〔公開〕〔規約〕

〔公開〕マニュアル 46.4.1 節：低バイアスの極限で \\( G = (e^2/h)\\,T(\mu\_{\mathrm L}) \\) を、\\( G\_0 = e^2/h \\) 単位でスピンごとに出力します。

〔規約〕実装は \\( \bar\mu = (\mu\_{\mathrm L} + \mu\_{\mathrm R})/2 \\) での透過率を、透過率のエネルギー格子の隣り合う 2 点の線形補間で求めます
（バイアス 0 ではマニュアルの \\( \mu\_{\mathrm L} \\) と同じ）。k 点ごとの値と平均を `System.Name.conductance` に出します。

### 9.4 バイアス効果の補間〔公開〕〔規約〕

〔公開〕技術ノート 13 Eq. 54–56、マニュアル 46.7 節：2 つのバイアスで SCF を済ませ、\\( H\_{\mathrm C} \\)、\\( H\_{\mathrm R} \\)、\\( \mu\_{\mathrm R} \\) を
\\( \lambda \\) と \\( 1 - \lambda \\) で線形補間して透過率・電流を求めます（`NEGF.tran.interpolate`、`NEGF.tran.interpolate.file1`、`file2`、`NEGF.tran.interpolate.coes`）。

〔規約〕実装は C の \\( H \\)、左右両方の電極の \\( H \\)、左右の化学ポテンシャル、C の化学ポテンシャルのすべてを同じ係数で補間します
（\\( V\_{\mathrm b} \\) を左右に振り分けるため、左電極も変わる）。

### 9.5 力・応力・全エネルギー〔規約〕

NEGF の段階ではエネルギー重み付き密度行列を計算せず、バンドエネルギーの項も 0 とします。マニュアル 46.11 節の自動試験は全エネルギーと力を比較していますが、
NEGF で得られる力・全エネルギーの意味（開放系での変分性）は公開資料に記載がありません。この章では対象外とし、〔未確定〕に挙げます。

## 10. 固有チャネル〔公開〕〔規約〕

〔公開〕技術ノート 20（Eq. 8–33）は、Paulsson–Brandbyge の方法（PRB **76**, 115117 (2007)）を非直交基底に適用します。
Löwdin 直交化 \\( S\_{\mathrm C} = S^{1/2}S^{1/2\dagger} \\) のもとで

$$
\bar A\_{\mathrm L} = S^{1/2\dagger}\\, G\Gamma\_{\mathrm L}G^\dagger\\, S^{1/2},\qquad
\bar\Gamma\_{\mathrm R} = S^{-1/2}\\,\Gamma\_{\mathrm R}\\,S^{-1/2\dagger},\qquad
\hat T = \bar A\_{\mathrm L}\bar\Gamma\_{\mathrm R}
$$

とし、\\( \bar A\_{\mathrm L} = \sum\_n a\_n\\, u\_n u\_n^\dagger \\) から \\( \bar A\_{\mathrm L}^{1/2} = (\sqrt{a\_1}u\_1, \sqrt{a\_2}u\_2, \dots) \\) を作り、
Hermite 行列 \\( \bar A\_{\mathrm L}^{1/2\dagger}\bar\Gamma\_{\mathrm R}\bar A\_{\mathrm L}^{1/2} \\) を対角化して透過固有値 \\( t\_n \\) と直交する \\( \tilde t\_n \\) を得ます。
固有チャネルは \\( t\_n' = S^{-1/2\dagger}\bar A\_{\mathrm L}^{1/2}\tilde t\_n \\) を元の基底での係数として実空間に描きます。
\\( \sum\_n t\_n = T \\)、\\( 0 \le t\_n \le 1 \\) です（`tools/check_negf.py` で確認）。

〔公開〕マニュアル 46.4.3 節：`NEGF.tran.Channel`（既定 `on`）、k 点は `NEGF.Channel.Nkpoint` と `<NEGF.Channel.kpoint`（面内の分数座標）、
エネルギーは `NEGF.Channel.Nenergy` と `<NEGF.Channel.energy`（\\( \mu\_{\mathrm L} \\) 基準の eV）、実空間に描く数は `NEGF.Channel.Num`（共線 5、非共線 10）。
透過固有値は大きい順に出力します。

〔規約〕

- \\( S^{1/2} \\) は \\( S\_{\mathrm C} \\) の固有分解から作り、固有値 \\( 10^{-6} \\) 未満の方向を除きます（標準出力の `N_{ort} / N_{nonort}` が残した次元と元の次元）。
  除いた方向のチャネルは固有値 0 として埋めます。
- 虚部は透過率と同じ \\( \eta \\)（`NEGF.tran.energyrange` の 3 番目）を使います。
- 各チャネルの全体の位相は、係数の和（虚部の符号で向きをそろえたもの）が実数になるように選びます。位相は物理量に影響しない慣習で、照合は \\( \lvert t'\_n\rvert \\) か密度 \\( \lvert T\_n(\mathbf r)\rvert^2 \\) で行います〔設計〕。

## 11. 電流密度〔公開〕〔規約〕

〔公開〕技術ノート 20（Eq. 34–58）は、電流密度を局所項と非局所項に分けます。
局所項は \\( D\_{\mathrm R} = G\Gamma\_{\mathrm R}G^\dagger \\)（非直交基底の行列）を使って

$$
\mathbf J\_{\mathrm{Loc}}(\mathbf r) = \int\frac{dE}{2\pi}\left[f(E - \mu\_{\mathrm L}) - f(E - \mu\_{\mathrm R})\right] i \sum\_{\mathbf q}\sum\_{ij} D^{\mathbf q}\_{ij}
\left\\{2i\mathbf k\\,\chi\_j^\*\chi\_i + \chi\_j^\*\nabla\chi\_i - \chi\_i\nabla\chi\_j^\*\right\\}
$$

（技術ノート 20 Eq. 36–40, 64）で、非局所擬ポテンシャルによる非局所項は、電荷保存から決まる Poisson 方程式 \\( \nabla^2\varphi = -\rho\_{\mathrm{NL}} \\)、
\\( \mathbf J\_{\mathrm{NL}} = -\nabla\varphi \\) の解として求めます（Eq. 41–47）。境界では Landauer–Büttiker の電流から Neumann 条件を与えます（Eq. 48–55）。
電荷電流は \\( \mathbf J\_\uparrow + \mathbf J\_\downarrow \\)、スピン電流は \\( \mathbf J\_\uparrow - \mathbf J\_\downarrow \\) です（Eq. 57–58）。
非共線ではスピンの 2×2 成分ごとに定義します（Eq. 63–72）。

〔公開〕マニュアル 46.4.2 節：`NEGF.tran.CurrentDensity`（既定 `on`）、非共線のスピン非対角成分は `NEGF.OffDiagonalCurrent`（既定 `off`）。

〔規約〕

- エネルギー積分は 9.2 節の電流と同じ格子と虚部で行います。\\( \lvert\mu\_{\mathrm L} - \mu\_{\mathrm R}\rvert < 10^{-6} \\) Hartree では計算しません。
- \\( S^{-1} \\) は \\( S\_{\mathrm C} \\) の固有分解で、固有値 \\( 10^{-6} \\) 未満の方向を除いた擬似逆行列とします。
- 標準出力の `Left(ideal)`／`Right(ideal)` は Landauer–Büttiker の境界電流の面積分、`Left(truncated)`／`Right(truncated)` は
  格子上の電流密度を境界面で積分した値です。両者の差が電流密度の離散化と基底の打ち切りの誤差の目安になります。
- 各原子の Voronoi 多面体での平均を `_atom` のファイルに出します。

詳細な離散化（勾配の評価、Neumann 条件の実装）は本版では書いていません（〔未確定〕）。

## 12. 非共線 NEGF〔公開〕〔規約〕

〔公開〕マニュアル 46.9 節：非共線 DFT と組み合わせられ、スピン軌道相互作用、DFT+U、スピン・軌道磁気モーメントの向きの拘束と併用できます。
手順は共線と同じです。

〔規約〕

- C と電極の行列は \\( 2N \\) 次元で、並びはスピンが外側（[k 点](kpoints-bloch.md) の非共線スピノルの並び）、重なりは \\( I\_2\otimes S \\) です。
  \\( \downarrow\uparrow \\) ブロックは \\( \uparrow\downarrow \\) の Hermite 共役で埋めます。
- k 点の時間反転による削減はしません（SCF・Step 3 とも）。
- バイアスは \\( \uparrow\uparrow \\)、\\( \downarrow\downarrow \\) にだけ \\( V\_sS \\) を加えます（6.1 節）。
- 透過率のファイルはスピンを分けず、全透過率 1 つの欄（実部・虚部）です。`NEGF.Channel.Num` の既定は 10 です。
- 境界付近の密度の補い（5.5 節）では、密度行列の \\( \uparrow\downarrow \\) 成分の虚部の符号を反転して足します。

## 13. その他の機能〔公開〕〔規約〕

- **Step 3 だけの再計算**（マニュアル 46.5 節）：`NEGF.tran.SCF.skip on` で `System.Name.tranb` を読んで Step 3 だけを行います。
  Step 2 だけで止めるには `NEGF.tran.SCF.skip off`、`NEGF.tran.Analysis off`、`NEGF.tran.Channel off` とします。
  〔規約〕`NEGF.tran.Channel`、`NEGF.tran.CurrentDensity`、`NEGF.OffDiagonalCurrent` のどれかが `on` なら `NEGF.tran.Analysis` は自動的に `on` になります。
  `NEGF.tran.SCF.skip` には `Periodic` という値もあります（次項の三次元周期系用）。
- **バイアス 0 の周期系**（マニュアル 46.6 節）：Step 1 に `NEGF.Output.for.TranMain on` を加えると、a 方向にも周期的な系の
  通常のバンド計算の \\( H \\) から `*.tranb` を作り、Step 2 を省いて透過率を計算できます。
- **状態密度**（マニュアル 46.3 節 E）：`Dos.fileout on` のとき、`NEGF.Dos.energyrange`（下限・上限・虚部、既定 \\( -10, 10, 5\times 10^{-3} \\) eV）、
  `NEGF.Dos.energy.div`（既定 200）、`NEGF.Dos.Kgrid`（既定 `1 1`）で \\( -\pi^{-1}\mathrm{Im}\\,G \\) から状態密度を作り、
  `System.Name.Dos.val`、`System.Name.Dos.vec` に出します。エネルギーの原点は \\( \mu\_{\mathrm L} \\) です。
- **並列化**（マニュアル 46.8 節）：エネルギー点・k 点・スピンの組を MPI で分け、行列演算を OpenMP で並列化します。結果には影響しません。
- **併用できないもの**〔規約〕：基底の最適化（`orbitalOpt`）とは同時に使えず、停止します。

## 14. 出力ファイル〔公開〕〔規約〕

| ファイル | 内容 | 出典 |
|---|---|---|
| `NEGF.filename.hks` で指定（Step 1） | 電極データ（2.1 節） | マニュアル 46.2 |
| `System.Name.tranb`（Step 2） | Step 3 に必要な C と電極の行列・化学ポテンシャル・領域の印など | マニュアル 46.1, 46.5 |
| `System.Name.paranegf` | 電極の化学ポテンシャル（固有値とバイアス後）、正味のバイアス、非平衡積分の区間・刻み・点数。`System.Name.out` に統合される〔規約〕 | 〔規約〕 |
| `System.Name.tran#_%` | k 点 \\( (\\#, \\%) \\) ごとの透過率。列は番号、\\( E - \mu\_{\mathrm L} \\) と虚部（Hartree）、同（eV）、\\( T\_\uparrow \\) の実部・虚部、\\( T\_\downarrow \\) の実部・虚部 | マニュアル 46.4.1、〔規約〕 |
| `System.Name.current` | k 点ごとの電流と平均（A、スピンごと） | マニュアル 46.4.1 |
| `System.Name.conductance` | k 点ごとのコンダクタンスと平均（\\( G\_0 = e^2/h \\) 単位、スピンごと） | マニュアル 46.4.1 |
| `System.Name.traneval#k_#E_#s`、`tranevec#k_#E_#s` | 透過固有値、固有チャネルの LCAO 係数 | マニュアル 46.4.3 |
| `System.Name.tranec#k_#E_#s_#c_r.cube`／`_i.cube` | 固有チャネルの実部・虚部（Gaussian cube） | マニュアル 46.4.3 |
| `System.Name.curden.xsf`、`scurden.xsf` | 電荷・スピン電流密度（格子上のベクトル場） | マニュアル 46.4.2 |
| `System.Name.curden1.cube`、`scurden1.cube` | 同、a 成分 | マニュアル 46.4.2 |
| `System.Name.curden_atom.xsf`、`scurden_atom.xsf` | 同、原子の Voronoi 多面体での平均 | マニュアル 46.4.2 |
| `System.Name.odcurden*_r`／`_i` | 非共線のスピン非対角成分（実部・虚部） | マニュアル 46.4.2 |

\\( \\# \\)、\\( \\% \\) は面内の k 点の番号（0 始まり）です。

## 15. 入力キーワードと既定値（まとめ）

〔公開〕はマニュアル 46 節に記載のもの、〔規約〕は記載がないか既定値が異なるものです。キーワードは大文字・小文字を区別しません〔規約〕。

| キーワード | 既定値 | 意味 | 区分 |
|---|---|---|---|
| `scf.EigenvalueSolver` | — | `NEGF` で Step 2 | 公開 |
| `NEGF.output_hks`、`NEGF.filename.hks` | `off`、`NEGF.hks` | Step 1 で電極データを書く | 公開 |
| `NEGF.filename.hks.l`、`.r` | `NEGF.hks.l`、`NEGF.hks.r` | 左右の電極データ | 公開 |
| `LeftLeadAtoms.*`、`RightLeadAtoms.*` | — | \\( \mathrm L\_0 \\)、\\( \mathrm R\_0 \\) の原子 | 公開 |
| `NEGF.Num.Poles` | 150 | 極の数 \\( N\_{\mathrm p} \\) | 公開 |
| `NEGF.scf.Kgrid` | `1 1` | SCF の面内 k 点 | 公開 |
| `NEGF.scf.Iter.Band` | 10 | 最初に Band で回す反復数（マニュアルは 6） | 規約 |
| `NEGF.ChemP.Band` | `fix` | Band の段階の化学ポテンシャル | 規約 |
| `NEGF.Band2NEGF.Hartree` | `off` | 切り替え時の Hartree の混合（7.3 節） | 規約 |
| `NEGF.bias.apply` | `on` | バイアス処理の有無 | 規約 |
| `NEGF.bias.voltage` | 0.0 eV | \\( V\_{\mathrm b} \\) | 公開 |
| `NEGF.bias.neq.im.energy` | 0.01 eV | \\( \eta\_{\mathrm{neq}} \\)（負なら停止） | 公開 |
| `NEGF.bias.neq.energy.step` | 0.02 eV | \\( \Delta E \\)（負なら停止） | 公開 |
| `NEGF.bias.neq.cutoff` | \\( 10^{-8} \\) | 非平衡の窓の打ち切り | 規約 |
| `NEGF.NonSCF.bias` | `off` | 6.3 節 | 規約 |
| `NEGF.gate.voltage` | 0.0 eV | \\( V\_{\mathrm g}^{(0)} \\) | 公開 |
| `NEGF.Poisson.Solver` | `FD` | `FD`／`FFT`（／`FDG`） | 公開（`FDG` は規約） |
| `NEGF.Poisson_Gparazero.scaling` | 1.0 | 7.1 節 | 規約 |
| `NEGF.Mixing.H.Boundary` | \\( -1 \\)（tanh の重み） | 2.3 節 | 規約 |
| `NEGF.Integration` | `CF` | `CF`／`OLD` | 規約 |
| `NEGF.Surfgreen.iterationmax` | 600 | 繰り込みの最大回数 | 規約 |
| `NEGF.Surfgreen.convergeeps` | \\( 10^{-12} \\) | 繰り込みの収束判定 | 規約 |
| `NEGF.tran.Analysis` | `on` | Step 3 の透過率・電流 | 公開 |
| `NEGF.tran.SCF.skip` | `off` | Step 3 だけ行う | 公開 |
| `NEGF.tran.energyrange` | \\( -10, 10, 10^{-3} \\) eV | 透過率の範囲と虚部 | 公開 |
| `NEGF.tran.energydiv` | 200 | 透過率の点数 | 公開 |
| `NEGF.tran.Kgrid` | `NEGF.scf.Kgrid` の値 | Step 3 の面内 k 点 | 公開（既定は規約） |
| `NEGF.current.energy.step` | 0.01 eV | 電流の積分の刻み | 規約 |
| `NEGF.current.im.energy` | `NEGF.tran.energyrange` の虚部 | 電流の積分の虚部 | 規約 |
| `NEGF.current.cutoff` | \\( 10^{-8} \\) | 電流の窓の打ち切り | 規約 |
| `NEGF.tran.Channel` | `on` | 固有チャネル | 公開 |
| `NEGF.Channel.Nkpoint`、`<NEGF.Channel.kpoint` | 1、\\( (0, 0) \\) | 固有チャネルの k 点 | 公開 |
| `NEGF.Channel.Nenergy`、`<NEGF.Channel.energy` | 1、0 eV | 固有チャネルのエネルギー | 公開 |
| `NEGF.Channel.Num` | 5（非共線 10） | 実空間に描くチャネル数 | 公開 |
| `NEGF.tran.CurrentDensity` | `on` | 電流密度 | 公開 |
| `NEGF.OffDiagonalCurrent` | `off` | 非共線のスピン非対角電流 | 公開 |
| `NEGF.tran.interpolate`（`.file1`、`.file2`、`.coes`） | `off`（係数 `1.0 0.0`） | バイアスの補間 | 公開 |
| `NEGF.Output.for.TranMain` | `off` | 三次元周期系から tranb を作る | 公開 |
| `NEGF.Dos.energyrange`、`.energy.div`、`.Kgrid` | \\( -10, 10, 5\times 10^{-3} \\) eV、200、`1 1` | NEGF の状態密度 | 公開 |
| `NEGF.FFTE.Num.Buffer.Cells` | 1 | 未使用とみられる（〔未確定〕） | 規約 |

## 16. 端の場合と注意〔規約〕〔設計〕

- **バイアス 0**：非平衡の点は 0 個で、平衡部分だけになります（5.3 節）。電流は 0、電流密度は計算しません。
- **バイアスの符号**：\\( V\_{\mathrm b} < 0 \\) でも \\( \mu\_<, \mu\_> \\) を入れ替えて同じ式を使います。参照電極は常に化学ポテンシャルの低い方です。
- **表面 Green 関数の非収束**：600 回で収束しない場合も計算は続きます。バンド端のごく近くで \\( \eta \\) が小さいと起こり得ます。
  再実装では非収束を警告として数え、報告します〔設計〕。
- **\\( \mathrm C\_0 \\) が電極と同じ**：完全な結晶では、透過率は面内の各 k 点でチャネル数（整数）になります（テスト参照）。
- **区間の刻みが粗い**：\\( E\_{\max} - E\_{\min} < \Delta E \\) なら警告を出し、非平衡部分の点が 0 個になります。
- **格子の a 方向の間隔**：C の \\( N\_1 \\) は電極の格子間隔と独立に決まるので、FD の境界値（電極の格子面の値）と C の格子面の間隔が厳密には一致しません。

## 実装への要求〔設計〕

- **モジュールの分割**：(1) 電極データ（hks に相当する、実空間の \\( H, S, \rho \\)、格子上の \\( \delta V\_{\mathrm H} \\)、密度、化学ポテンシャル）の型、
  (2) 拡張中心領域の組み立て（原子の対応、領域の印、境界の混合）、(3) k 点ごとのブロック（\\( H\_{00}, H\_{\mathrm{out}}, H\_{\mathrm Cs}, H\_{\mathrm C} \\)）の生成、
  (4) 表面 Green 関数、(5) 自己エネルギーと \\( G \\)、(6) 積分点と重みの生成（極展開・モーメント・非平衡の矩形則）、
  (7) 密度行列の組み立て（k 点の和）、(8) 境界条件付き Poisson ソルバー、(9) Step 3 の解析（透過率・電流・固有チャネル・電流密度）、に分けます。
- **積分点の表現**：積分点を「複素エネルギー \\( z \\)、重み \\( w \\)、種類（極・モーメント・非平衡）」の列として生成し、
  密度行列は \\( \sum w\\,X(z) \\)（\\( X \\) は \\( G + G^\dagger \\)、\\( G \\)、\\( G\Gamma G^\dagger \\) のいずれか）として足し込みます。
  点ごとに独立なので、点・k 点・スピンの組で並列化できます。
- **極展開**：\\( N\_{\mathrm p} \\) から \\( z\_p, R\_p \\) を返す関数を独立させ、5.1 節の式で検査できるようにします。
- **バイアスの振り分け**：既存実装との照合のため、\\( (V\_{\mathrm L}, V\_{\mathrm R}) = (-V\_{\mathrm b}/2, +V\_{\mathrm b}/2) \\) を既定とし、引数で変えられるようにします。
- **表面 Green 関数**：繰り込み法を既定とし、収束判定値・最大回数・実際の回数・最終の差を返します。Dyson 残差を検査用に計算できるようにします。
- **自己エネルギーの再利用**：\\( g\_s(z) \\) は電極と \\( z \\)、\\( \mathbf q \\) だけで決まり C の SCF に依らないので、SCF の反復をまたいでキャッシュできます
  （既存実装は毎回計算し直しています。結果は変わりません）。
- **0 次モーメント**：\\( iR\\,G(iR) \\) を自己エネルギー込みで計算します。\\( R \\) は引数（既定 \\( 10^{12} \\) Hartree）。
- **k 点**：面内の `Regular` 格子と時間反転の削減は [k 点](kpoints-bloch.md) の関数を二次元で使います。
- **単位**：入力の eV は Hartree に直して保持し、出力の電流は \\( eE\_{\mathrm h}/\hbar \\) を定数として明示します（既存実装の値と CODATA の値を選べるようにする）。

## テスト〔設計〕

`tools/check_negf.py` は、2 軌道の非直交な電極胞を持つ一次元模型で以下の多くを確認しています。

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 極展開（\\( N\_{\mathrm p} = 150 \\)）の \\( \lvert x\rvert \le 60 \\) での最大誤差 | \\( 10^{-12} \\) |
| 解析 | 最初の 20 個の \\( z\_p \\) が \\( (2p-1)\pi \\)、\\( R\_p = -1 \\) | \\( 10^{-10} \\)（相対） |
| 解析 | 1 軌道の鎖（\\( \varepsilon\_0, t \\)）の表面 Green 関数 \\( g = [w - \sqrt{w - 2t}\sqrt{w + 2t}]/(2t^2) \\)、\\( w = z - \varepsilon\_0 \\) | \\( 10^{-7} \\) |
| 不変量 | 繰り込みの結果が表面 Green 関数の Dyson 方程式を満たす | \\( 10^{-9} \\) |
| 不変量 | 完全な鎖（\\( \mathrm C\_0 \\) = 電極）で \\( T(E) \\) = 電極のバンドから数えたチャネル数（\\( \eta = 10^{-8} \\)） | \\( 10^{-5} \\) |
| 不変量 | \\( T\_{\mathrm{LR}} = \mathrm{Tr}[\Gamma\_{\mathrm L}G\Gamma\_{\mathrm R}G^\dagger] = T\_{\mathrm{RL}} \\) | \\( 10^{-10} \\) |
| 不変量 | 透過固有値の和 = \\( T \\)、\\( 0 \le t\_n \le 1 \\) | \\( 10^{-10} \\) |
| 不変量 | 孤立系で、極展開の密度行列 = \\( \sum\_n c\_n f(\varepsilon\_n) c\_n^\dagger \\) | \\( 10^{-12} \\) |
| 不変量 | 電極付きの系で、\\( N\_{\mathrm p} = 150 \\) と 300 の密度行列が一致、密度行列が Hermite | \\( 10^{-11} \\) |
| 照合（独立計算） | 極展開の密度行列と実軸積分（\\( \eta = 2\times 10^{-3} \\)）の差 | \\( 5\times 10^{-3} \\)（\\( O(\eta) \\)） |
| 不変量 | 非平衡の密度行列が参照電極（\\( \mu\_< \\) と \\( \mu\_> \\)）によらない | \\( 2\times 10^{-3} \\)（\\( \eta\_{\mathrm{neq}} \\)、刻みに依存） |
| 解析 | \\( eE\_{\mathrm h}/\hbar = 6.6236\times 10^{-3} \\) A | \\( 10^{-6} \\)（相対） |
| 不変量 | 時間反転の対の k 点で \\( T(E; \mathbf q) = T(E; -\mathbf q) \\)（共線） | \\( 10^{-10} \\) |
| 不変量 | 左右を入れ替えた系（鏡映）で \\( T \\) が不変、\\( V\_{\mathrm b} \to -V\_{\mathrm b} \\) で電流の符号が反転 | \\( 10^{-8} \\) |
| 不変量 | 電極の全ブロックに \\( VS \\) を加え \\( \mu \\) も \\( V \\) ずらすと、\\( T(E + V) \\) が元の \\( T(E) \\) に一致（C も同じだけずらした場合） | \\( 10^{-10} \\) |
| 照合 | マニュアルの例題（炭素鎖 Lead-Chain.dat／NEGF-Chain.dat、グラフェン、8-ZGNR の有限バイアス 0.3 V、非共線 8-ZGNR、非共線の金の鎖）の Step 2 の `NormRD` の収束後の Mulliken 電荷 | \\( 10^{-6} \\) |
| 照合 | 同じ例題の透過率曲線（`System.Name.tran0_0` の 6・8 列） | \\( 10^{-6} \\)（既存実装の \\( \eta \\) を合わせる） |
| 照合 | 同じ例題の電流、コンダクタンス、透過固有値の和 | \\( 10^{-7} \\)（マニュアル 46.11 節の自動試験は小数点以下 7 桁で比較） |
| 照合 | 既存実装の `System.Name.paranegf` の非平衡積分の区間・点数を、5.3 節の規則で再現 | 点数は厳密、区間は \\( \Delta E/10 \\) |

照合用の fixture は、既存実装で各例題の Step 1・Step 2 を実行し、hks（電極データ）、tranb、標準出力、Step 3 の出力ファイルを保存して作ります。
並列数で結果がわずかに変わり得るので、MPI の数を記録します。

## 〔未確定〕

- `NEGF.bias.apply off` の経路：既存実装はこの経路で極の重みを \\( -2R\_p/\beta \\)（既定の経路の 2 倍）にしており、モーメント項は同じ \\( \tfrac12 iR\\,G(iR) \\) のままです。
  5.1 節の式と一致しないように見えます。意図的なもの（別の箇所の因子との組み合わせ）か不具合かを、実行して確認する必要があります。
- NEGF での力・応力・全エネルギーの意味と計算方法（エネルギー重み付き密度行列を計算していない。技術ノート 13 の付録は非平衡のエネルギー密度行列の式を与えている）。
- NEGF 用の電荷ミキシング（bc 面の Fourier 成分と a 方向の格子面の表現での RMM-DIIS・Kerker）の詳細。
- 電流密度の離散化（勾配、非局所項の Poisson 方程式と Neumann 境界条件）の実装の詳細と、`ideal`／`truncated` の定義の厳密な内容。
- `FDG` 方式と 2 値の `NEGF.gate.voltage`、`NEGF.FFTE.Num.Buffer.Cells` の意味。
- C の \\( N\_1 \\) と電極の a 方向の格子間隔が一致しない場合の FD 境界条件の誤差の評価。
- 非共線 NEGF での電流・コンダクタンスのファイル形式の詳細（スピンの欄の扱い）、`NEGF.OffDiagonalCurrent` の出力の規約。
- `System.Name.tranb` と hks ファイルのバイト配置（既存実装の出力を読む互換層を作るかどうか）。
- 参照したのは GPU 版の OpenMX 4.0 で、本家の CPU 版で既定値（特に `NEGF.scf.Iter.Band`）が同じかは確認していません。
