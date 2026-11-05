# 構造最適化・分子動力学・NEB

> 状態：初版（レビュー待ち）。`MD.Type` のすべての選択肢（構造最適化、格子の最適化、分子動力学、NEB）について、
> 公開資料（マニュアル 17–20・51 節、技術ノート 7「Geometry Optimization」、NEB の技術ノート、標準的な論文）の式を
> このノートの記法で書き直し、公開資料にない規約・既定値を既存実装から確認して加えたものです。
> 数値で確かめた部分は `tools/check_geometry_md.py` にあります。

## この章の範囲

- SCF が収束して力（と応力）が求まった後に、原子座標・格子ベクトル・速度を次のステップへ進める部分を扱います。
- 力と応力の定義・符号は [全エネルギーと力](energy-force.md)、次のステップの SCF の初期密度（外挿）は
  [密度行列と電荷ミキシング](density-mixing.md) の 6・7 節、格子の原点は [実空間グリッド](grid-poisson-xc.md) の 1.4 節にあります。
  この章ではそれらへの入口だけを書きます。
- `MD.Type` の値ごとに、何を変数にし、どの式で更新し、何を出力するかを定めます。
  再実装の目標は「同じ入力で、同じ構造に収束する」「同じ統計集団の軌跡を生成する」ことで、
  既存実装とステップごとに一致させることは、明記した場合（照合用の互換モード）に限ります。

## 記号〔設計〕

[記法](../design/notation.md) の記号に加えて、次を使います。内部単位は原子単位（bohr、Hartree、電子質量 \\( m\_e \\)、
時間の単位 \\( \hbar/E\_h \\)）です。

| 記号 | 意味 |
|---|---|
| \\( N \\) | 原子数（単位胞内） |
| \\( \mathbf R\_I \\)、\\( \mathbf x \in \mathbb R^{3N} \\) | 原子 \\( I \\) の位置と、全原子の座標を並べたベクトル |
| \\( U \\) | 全エネルギー `Utot`（[全エネルギーと力](energy-force.md)） |
| \\( \mathbf g = \partial U/\partial\mathbf x \\)、\\( \mathbf F = -\mathbf g \\) | 勾配と力（Hartree/bohr） |
| \\( m\_I \\) | 原子 \\( I \\) の質量（内部では \\( m\_e \\) 単位） |
| \\( \mathbf v\_I \\)、\\( K = \sum\_I \tfrac12 m\_I \lvert\mathbf v\_I\rvert^2 \\) | 速度と原子核の運動エネルギー `Ukc` |
| \\( \Delta t \\)、\\( n \\) | 時間刻みと、ステップ番号（1 始まり。再開時は通算） |
| \\( A \\) | 格子ベクトルを **行** とする 3×3 行列（\\( A\_{ij} \\) は \\( \mathbf a\_i \\) の \\( j \\) 成分） |
| \\( \Omega = \det A \\) | 単位胞の体積 |
| \\( \boldsymbol\sigma\_i \\)、\\( S \\) | \\( \boldsymbol\sigma\_1 = \mathbf a\_2\times\mathbf a\_3 \\) など（巡回）。これを行とする行列 \\( S = \Omega A^{-\mathsf T} \\) |
| \\( \mathbf s\_I \\) | 分数座標 \\( \mathbf R\_I = \mathbf O + A^{\mathsf T}\mathbf s\_I \\)（\\( \mathbf O \\) は格子の原点） |
| \\( \Sigma\_{\gamma\eta} \\) | [全エネルギーと力](energy-force.md) の応力（\\( \Sigma/\Omega \\) が出力される応力テンソル。引っ張られた状態で対角が正） |
| \\( k\_B \\)、\\( T \\)、\\( T\_0 \\) | Boltzmann 定数、運動学的温度、目標温度 |
| \\( \zeta \\)、\\( \eta \\)、\\( Q \\) | Nosé–Hoover の摩擦係数、その時間積分（熱浴の座標）、熱浴の質量 |
| \\( W \\)、\\( P \\) | バロスタットの質量と外圧 |

## 1. 全体の流れ〔公開〕〔規約〕

〔公開〕マニュアル 18・20 節：構造最適化と MD は同じ入力ファイルで `MD.Type` を選び、ステップごとに
SCF → 力 → 座標の更新 を繰り返します。

〔規約〕1 ステップの処理は次の順です。

1. 現在の座標（と格子）で SCF を収束させ、`Utot`、力、必要なら応力を求める。
2. `.md` ファイルに現在の座標・力（と、その時点で保持している速度）を 1 ブロック追記する（10 節）。
3. 分子動力学では、力の総和を 0 にする補正をかける（下記）。
4. `MD.Type` に応じて座標・格子・速度を更新し、`.ene` ファイルに 1 行追記する。
5. 分子動力学（`NVE`、`NVT_*`、`NPT_*`）に限り、単位胞の外に出た原子を中に戻す（7 節）。構造最適化では戻しません。
6. 現在の状態から再開用の入力ファイル `System.Name.dat#` を書き出す（10 節）。

- ステップ番号は通算の \\( n = n\_{\mathrm{local}} + \\) `MD.Current.Iter` です（`MD.Current.Iter` は既定 0、再開ファイルが書き込む）。
  繰り返しは、収束の判定が立つか \\( n \\) が `MD.maxIter`（既定 1）に達したところで終わります。
  `MD.maxIter` は再開前のステップを含めた **通算** の上限です。
- \\( n = \\) `MD.maxIter` のステップでは、収束していなくても座標を更新しません（最後の SCF の構造と力が対応したまま終わる）。
- 分子動力学（`NVE`、`NVT_VS`、`NVT_VS2`、`NVT_VS4`、`NVT_NH`、`NVT_Langevin`）では、更新の前に、
  固定されていない成分について、デカルト成分ごとに勾配の **単純平均**（質量で重みを付けない）を引きます：
  \\( g\_{I\alpha} \leftarrow g\_{I\alpha} - \frac{1}{N\_\alpha}\sum\_{J\\,\text{free}} g\_{J\alpha} \\)（\\( N\_\alpha \\) は \\( \alpha \\) 成分が自由な原子の数）。
  NPT と構造最適化ではこの補正をしません（[全エネルギーと力](energy-force.md) の「符号と単位」と同じ内容）。
- 2 ステップ目以降の SCF は、前のステップの密度から再開し、密度を外挿します（[密度行列と電荷ミキシング](density-mixing.md) の 6・7 節）。
  格子の原点 \\( \mathbf O \\) は最初のステップで一度決めたら変えません（[実空間グリッド](grid-poisson-xc.md) 1.4 節）。

## 2. 単位と定数〔規約〕

| 量 | 入力の単位 | 内部単位への換算 |
|---|---|---|
| 時間刻み `MD.TimeStep` | fs（既定 0.5、負なら入力エラー） | \\( 1\ \mathrm{fs} = 41.3411 \\) |
| 質量（`Definition.of.Atomic.Species` の 4 列目） | 統一原子質量単位 | \\( 1\ \mathrm u = 1.660538921/0.000910938291\ m\_e \approx 1822.888\ m\_e \\) |
| 速度（`MD.Init.Velocity`、`.md` の出力） | m/s | \\( 1\ \mathrm{m/s} = 0.4571028\times 10^{-6} \\) |
| 圧力（`MD.applied.pressure`、`MD.CellPressureControl`） | GPa | \\( 1\ \mathrm{GPa} = 0.3398827\times 10^{-4}\ \mathrm{Hartree/bohr^3} \\) |
| 温度 | K | \\( k\_B = 8.617251324\times10^{-5}\ \mathrm{eV/K} \\)、\\( 1\ \mathrm{Hartree} = 27.2113845\ \mathrm{eV} \\) |
| 長さ | Å または bohr | \\( 1\ \mathrm{bohr} = 0.529177249\ \text{Å} \\) |

- GPa の換算係数は CODATA の値（\\( 3.39893\times10^{-5} \\)）と相対 \\( 3\times10^{-5} \\) ずれています。照合では既存実装の係数を使います。
- 質量の既定値は、各元素で最も存在比の大きい同位体の原子質量です（例：H 1.00782503207、C 12.00、O 15.99491461956、Si 27.9769265325、Fe 55.9349375）。
  〔公開〕マニュアル 20.9 節：`Definition.of.Atomic.Species` の 4 列目で上書きできます。
- 〔規約〕温度の定義は 2 種類あり、場所によって使い分けられています。
  - \\( T\_A = K/(\tfrac32 k\_B N) \\)：全原子数 \\( N \\) で割る。`.ene` の温度、熱浴の制御、NPT で使う。
  - \\( T\_B = 2K/(k\_B N\_{\mathrm{free}}) \\)：\\( N\_{\mathrm{free}} \\) は固定されていない成分の数。
    `NVE`、`NVT_VS`、`NVT_NH`、`NVT_VS2`、`NVT_Langevin` の `.md` の見出し行の温度に使う。
  - どちらの場合も \\( K \\) は固定されていない成分だけで和を取ります（`NVT_VS` 系は全成分。6.4 節）。

## 3. `MD.Type` の一覧〔公開〕〔規約〕

〔公開〕マニュアル 18–20・51 節。既定値は `NOMD` です〔規約〕。

| 値 | 種類 | 変数 | 方法 | 節 |
|---|---|---|---|---|
| `NOMD` | 1 回の SCF | — | 座標を動かさない。`.ene` に 1 行、`.out` に座標と力を書いて終わる | — |
| `Opt` | 構造最適化 | 原子座標 | 最急降下（可変の前因子） | 4.2 |
| `DIIS` | 〃 | 〃 | GDIIS ＋ 単位行列の Hessian | 4.4 |
| `BFGS` | 〃 | 〃 | GDIIS ＋ BFGS の逆 Hessian | 4.5 |
| `RF` | 〃 | 〃 | GDIIS ＋ BFGS の Hessian ＋ 有理関数（RF）ステップ | 4.6 |
| `EF` | 〃 | 〃 | GDIIS ＋ BFGS の Hessian ＋ 固有値の下限処理 | 4.7 |
| `OptC1`–`OptC4` | 格子の最適化 | 格子（分数座標は固定） | 最急降下 | 5.2 |
| `OptC5`–`OptC7` | 〃 | 格子と原子座標 | 最急降下 | 5.2 |
| `RFC5`–`RFC7` | 〃 | 格子と分数座標 | GDIIS ＋ BFGS ＋ RF | 5.3 |
| `EvsLC` | エネルギー対格子定数 | 格子 | 一定の刻みで格子を伸縮 | 5.5 |
| `DF` | デルタ因子用 | 体積 | 7 点の体積で計算 | 5.5 |
| `NVE` | 分子動力学 | 座標・速度 | 速度 Verlet | 6.2 |
| `NVT_VS`、`NVT_VS2`、`NVT_VS4` | 〃 | 〃 | 速度スケーリング（全体・原子ごと・原子群ごと） | 6.4 |
| `NVT_NH` | 〃 | 〃 | Nosé–Hoover | 6.5 |
| `NVT_Langevin` | 〃 | 〃 | Langevin 熱浴 | 6.6 |
| `NPT_VS_PR`、`NPT_VS_WV`、`NPT_NH_PR`、`NPT_NH_WV` | 〃 | 座標・速度・格子 | 速度スケーリング／Nosé–Hoover と、Parrinello–Rahman／Wentzcovitch | 6.7 |
| `NEB` | 反応経路 | 全像の座標 | nudged elastic band ＋ DIIS・BFGS | 8 |

〔規約〕その他の扱い：

- `RFC1`–`RFC4` は受け付けるが「対応していない」と表示して終了します。
- `Constraint_DIIS`、`Opt_LBFGS` は受け付けますが、更新を何もしません（`MD.maxIter` まで同じ構造で SCF を繰り返す）。
  再実装では入力エラーにしてよい〔設計〕。
- `MD.Opt.DIIS.History` が 19 を超えると終了します。

## 4. 構造最適化（原子座標）

### 4.1 収束判定と最大変位〔公開〕〔規約〕

〔公開〕マニュアル 18.2 節：`MD.Opt.criterion`（既定 0.0003 Hartree/bohr）を、原子に働く力の最大値と比べます。

〔規約〕「力の最大値」の取り方が方式によって違います。固定された成分（4.10 節）は除きます。

| 方式 | 比べる量 |
|---|---|
| `Opt`（最急降下） | 原子ごとのベクトルの大きさの最大 \\( \max\_I \lVert\mathbf g\_I\rVert \\) |
| `DIIS`、`BFGS`、`RF`、`EF` | 成分ごとの絶対値の最大 \\( \max\_{I\alpha}\lvert g\_{I\alpha}\rvert \\) |

- 判定は、そのステップの SCF の後の勾配で行い、満たせばそのステップでは座標を更新せずに終わります。
- 1 ステップの最大変位：準 Newton 系では、更新量の成分ごとの最大 \\( \max\_{I\alpha}\lvert\Delta x\_{I\alpha}\rvert \\) が
  \\( d\_{\max} = 0.5 \\) bohr を超えたら、更新ベクトル全体を \\( d\_{\max}/\max\lvert\Delta x\rvert \\) 倍に縮めます（向きは保つ）。
  最急降下では \\( s\max\_I\lVert\mathbf g\_I\rVert \\)（下記）が 0.5 bohr を超えたら同様に縮めます。
- 〔規約〕GDIIS の 1 次方程式が解けない（係数行列が特異）場合、既存実装は **収束したものとして** 終了します。
  再実装では、履歴を捨てて最急降下に戻すか、エラーとして報告します〔設計〕。

### 4.2 最急降下法（`Opt`）〔公開〕〔規約〕

〔公開〕マニュアル 18.1 節：前因子が収束の様子に応じて変わる最急降下法です。

〔規約〕更新は \\( \mathbf x\_{n+1} = \mathbf x\_n - c\\, s\_n\\, \mathbf g\_n \\) で、\\( c \\) は 4.1 節の最大変位の縮小率です。前因子 \\( s \\)（bohr²/Hartree）は次のとおりです。

- 基準値 \\( s\_0 = (2\ \mathrm{fs})^2/(1\ \mathrm u) = (82.6822)^2/1822.888 \approx 3.750 \\)。範囲を \\( [0.04 s\_0,\ 15 s\_0] \\) に制限します。
- 最初のステップ（と、SD に切り替わった最初のステップ、4.8 節）：\\( s = 1/(1.5\times 0.529177249) \approx 1.260 \\)
  （最大の力の原子がおよそ \\( \lVert\mathbf g\rVert \\) Å の 1/1.5 だけ動く値）。
- 以後、各ステップの `Utot` を \\( U\_n \\) として
  - \\( U\_n > U\_{n-1} \\) なら \\( s \leftarrow s/2 \\)。
  - そうでなく、\\( U\_{n-1} < U\_{n-2} \\) かつ \\( U\_n < U\_{n-1} \\) で、\\( n \bmod 4 = 1 \\) なら \\( s \leftarrow 2.5 s \\)。
  - 比較に使うエネルギーの履歴は 2 ステップ目から記録します（1 ステップ目の値は比較に入らない）。

### 4.3 GDIIS（共通部分）〔公開〕

〔公開〕技術ノート 7 の Eq. 6–10（Csaszar–Pulay の GDIIS、技術ノートでは RMM-DIIS と呼ぶ）。
直近 \\( p \\) ステップの座標 \\( \mathbf x^{(m)} \\) と勾配 \\( \mathbf g^{(m)} \\) から、
\\( \sum\_m c\_m = 1 \\) の条件の下で \\( \lVert\sum\_m c\_m\mathbf g^{(m)}\rVert^2 \\) を最小にする係数を求めます：

$$
\begin{pmatrix} B & \mathbf 1 \\\\ \mathbf 1^{\mathsf T} & 0 \end{pmatrix}
\begin{pmatrix} \mathbf c \\\\ \lambda/2 \end{pmatrix}
{}= \begin{pmatrix} \mathbf 0 \\\\ 1 \end{pmatrix},\qquad
B\_{mm'} = \mathbf g^{(m)}\cdot\mathbf g^{(m')}
$$

$$
\tilde{\mathbf x} = \sum\_m c\_m \mathbf x^{(m)},\qquad
\tilde{\mathbf g} = \sum\_m c\_m \mathbf g^{(m)},\qquad
\mathbf x\_{n+1} = \tilde{\mathbf x} + \Delta\mathbf x(\tilde{\mathbf g})
$$

〔規約〕

- \\( p = \min(k, \\) `MD.Opt.DIIS.History`\\( ) \\)。\\( k \\) は準 Newton 段階に入ってからのステップ数（4.8 節）。`MD.Opt.DIIS.History` の既定は 2。
  ただし、最急降下の段階の座標・勾配も履歴に入っています。
- \\( B \\) はその要素の絶対値の最大で割ってから解きます（\\( \mathbf c \\) は変わらない）。対称不定値の直接法で解きます。
- 固定された成分の勾配は 0 として履歴に入れます。
- 閉じた形は \\( \mathbf c = B^{-1}\mathbf 1/(\mathbf 1^{\mathsf T}B^{-1}\mathbf 1) \\) です（`tools/check_geometry_md.py` で確認）。

### 4.4 `DIIS`〔公開〕〔規約〕

〔公開〕技術ノート 7 Eq. 11：Hessian を単位行列で近似し \\( \mathbf x\_{n+1} = \tilde{\mathbf x} - \alpha\tilde{\mathbf g} \\)。

〔規約〕\\( \alpha \\)（bohr²/Hartree）は \\( r = \lVert\tilde{\mathbf g}\rVert/N \\)（3N 次元の Euclid ノルムを原子数で割る）で決めます。

| \\( r \\) | \\( > 10^{-2} \\) | \\( (10^{-3}, 10^{-2}] \\) | \\( (10^{-4}, 10^{-3}] \\) | \\( (10^{-5}, 10^{-4}] \\) | \\( (10^{-6}, 10^{-5}] \\) | \\( \le 10^{-6} \\) |
|---|---|---|---|---|---|---|
| \\( \alpha \\) | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 |

### 4.5 `BFGS`〔公開〕〔規約〕

〔公開〕技術ノート 7 Eq. 12–14：逆 Hessian \\( H^{-1} \\) を BFGS で更新し、\\( \mathbf x\_{n+1} = \tilde{\mathbf x} - H^{-1}\tilde{\mathbf g} \\)。
\\( \Delta\mathbf x = \mathbf x^{(n)} - \mathbf x^{(n-1)} \\)、\\( \Delta\mathbf g = \mathbf g^{(n)} - \mathbf g^{(n-1)} \\)（生の履歴の差。GDIIS の平均ではない）として

$$
H^{-1}\_{\mathrm{new}} = H^{-1} + \frac{\Delta\mathbf x\cdot\Delta\mathbf g + \Delta\mathbf g\cdot H^{-1}\Delta\mathbf g}{(\Delta\mathbf x\cdot\Delta\mathbf g)^2}\\, \Delta\mathbf x\\,\Delta\mathbf x^{\mathsf T}
{}- \frac{H^{-1}\Delta\mathbf g\\,\Delta\mathbf x^{\mathsf T} + \Delta\mathbf x\\,\Delta\mathbf g^{\mathsf T}H^{-1}}{\Delta\mathbf x\cdot\Delta\mathbf g}
$$

〔規約〕準 Newton 段階の最初のステップで \\( H^{-1} = I \\)（単位は bohr²/Hartree）とし、以後毎ステップ更新します。
\\( \Delta\mathbf x\cdot\Delta\mathbf g \le 0 \\) でも更新を止めません。`MD.Opt.Init.Hessian` は使いません（マニュアル 18.3 節も、モデル Hessian は BFGS・DIIS では使えないとしています）。

### 4.6 `RF`〔公開〕〔規約〕

〔公開〕技術ノート 7 Eq. 15–19（Banerjee ら、Baker の有理関数法）：Hessian \\( H \\) を BFGS で直接更新し

$$
H\_{\mathrm{new}} = H + \frac{\Delta\mathbf g\\,\Delta\mathbf g^{\mathsf T}}{\Delta\mathbf x\cdot\Delta\mathbf g} - \frac{H\Delta\mathbf x\\,(H\Delta\mathbf x)^{\mathsf T}}{\Delta\mathbf x\cdot H\Delta\mathbf x}
$$

拡大 Hessian の最小固有値 \\( \lambda \\) とその固有ベクトルから更新量を得ます：

$$
\begin{pmatrix} H & \tilde{\mathbf g} \\\\ \tilde{\mathbf g}^{\mathsf T} & 0\end{pmatrix}
\begin{pmatrix} \Delta\mathbf x \\\\ 1\end{pmatrix} = \lambda\begin{pmatrix} \Delta\mathbf x \\\\ 1\end{pmatrix},\qquad
\mathbf x\_{n+1} = \tilde{\mathbf x} + \Delta\mathbf x
$$

同値な形は \\( \Delta\mathbf x = -(H - \lambda I)^{-1}\tilde{\mathbf g} \\)、\\( \lambda = \tilde{\mathbf g}\cdot\Delta\mathbf x < \lambda\_{\min}(H) \\) です（`tools/check_geometry_md.py`）。
固有ベクトルは最後の成分で割って規格化します。

〔規約〕

- 初期 Hessian（4.9 節）は、準 Newton 段階の 1, 31, 61, … ステップ目（\\( k \bmod 30 = 1 \\)）に作り直します。
- BFGS の分母 \\( \Delta\mathbf x\cdot\Delta\mathbf g \\) と \\( \Delta\mathbf x\cdot H\Delta\mathbf x \\) は、\\( 10^{-10} \\) より小さければ \\( 10^{-10} \\) に置き換えて **更新を続けます**。
  曲率条件が破れた（\\( \Delta\mathbf x\cdot\Delta\mathbf g < 0 \\)）ステップでは、Hessian に大きな正の項 \\( 10^{10}\\,\Delta\mathbf g\Delta\mathbf g^{\mathsf T} \\) が入ります。
  マニュアル 18.2 節の「内積が負なら更新しない」という説明とは異なります。再実装では更新を飛ばす方を既定にします〔設計〕。

### 4.7 `EF`〔公開〕〔規約〕

〔公開〕技術ノート 7 Eq. 20–21（Baker の固有ベクトル追跡法を最小化に使う形）：BFGS で更新した \\( H \\) を対角化し、
固有値を正にそろえてから \\( H'^{-1} \\) を作り、\\( \mathbf x\_{n+1} = \tilde{\mathbf x} - f\\, H'^{-1}\tilde{\mathbf g} \\) とします。

〔規約〕

- 初期 Hessian（4.9 節）は準 Newton 段階の最初のステップでだけ作ります。
- BFGS の更新は \\( \Delta\mathbf x\cdot\Delta\mathbf g > 0 \\) のときだけ行います。
- 固有値の下限：\\( N \le 4 \\) なら 0.10、それ以外は 0.005 Hartree/bohr² より小さい固有値をその値に置き換えます
  （技術ノート 7 は Ver. 3.3 の値として 0.02 を挙げていますが、現在の値はこれです）。
- 前因子 \\( f \\)：準 Newton 段階の最初のステップで \\( f = 2.0 \\) とし、以後 `Utot` が前のステップより上がるたびに \\( f \leftarrow 0.95 f \\) とします。

### 4.8 最急降下との切り替え〔公開〕〔規約〕

〔公開〕マニュアル 18.2 節：`MD.Opt.StartDIIS` のステップまでは最急降下を使い、そこから準 Newton 法に切り替えます。
鞍点への収束を避けるには、最急降下で変曲点を越えてから切り替えるよう、この値を大きくします。

〔規約〕

| キーワード | 既定値（既存実装） | マニュアルの記載 | 意味 |
|---|---|---|---|
| `MD.Opt.DIIS.History` | 2 | 3 | GDIIS の履歴の長さ \\( p \\) の上限 |
| `MD.Opt.StartDIIS` | 1 | 5 | このステップから準 Newton 法 |
| `MD.Opt.EveryDIIS` | 200000 | 200 | 一方の段階がこの回数続いたら、もう一方に切り替える |
| `MD.Opt.criterion` | 0.0003 | 0.0003 | Hartree/bohr |
| `MD.maxIter` | 1 | 1 | 通算のステップ数の上限 |

- \\( n < \\) `MD.Opt.StartDIIS` の間は最急降下（4.2 節）。\\( n = \\) `MD.Opt.StartDIIS` で準 Newton 段階に入り、
  段階内のステップ番号 \\( k \\) を 1 から数えます（GDIIS の履歴の長さ、Hessian の作り直しはこの \\( k \\) で決まる）。
- 段階を切り替えたとき、最急降下の前因子は 4.2 節の初期値に戻ります。
- マニュアルの既定値の記載（表の 3 列目）は既存実装と一致しません。再実装は既存実装の値を既定にします〔設計〕。

### 4.9 初期 Hessian〔公開〕〔規約〕

〔公開〕マニュアル 18.3 節：`MD.Opt.Init.Hessian` は `Schlegel`（既定）または `iden`。
`Schlegel` は Schlegel のモデル（結合の伸縮・変角・ねじれ）のうち **結合の伸縮だけ** を使います。
モデルの出典は H. B. Schlegel, Theor. Chim. Acta **66**, 333 (1984)、J. M. Wittbrodt and H. B. Schlegel, J. Mol. Struct. (Theochem) **398–399**, 55 (1997)。

〔規約〕

- `iden`：\\( H = I \\)（Hartree/bohr²）。
- `Schlegel`：\\( H = 0.1\\, I \\) に、原子対の寄与を加えます。原子 \\( I \\) とその近接原子 \\( (J, \mathbf L) \\)
  （\\( \lVert\boldsymbol\tau\_{IJ\mathbf L}\rVert \le r\_{c,I} + r\_{c,J} \\)、[疎ブロック形式](../design/sparse-block-format.md) の原子対の集合。\\( (J,\mathbf L) \ne (I, 0) \\)）について、
  \\( r = \lVert\mathbf R\_I - \mathbf R\_J - \mathbf L\rVert \\)、\\( \mathbf u = (\mathbf R\_I - \mathbf R\_J - \mathbf L)/r \\)、

$$
k\_{IJ} = \frac{1.734}{(r - B\_{p(I)p(J)})^3}\ \ (\text{Hartree/bohr}^2),\qquad
H\_{I\alpha,I\beta} \mathrel{+}= k\_{IJ} u\_\alpha u\_\beta,\qquad
H\_{I\alpha,J\beta} \mathrel{+}= -k\_{IJ} u\_\alpha u\_\beta
$$

  （すべての順序対 \\( (I,J) \\) について加えるので、行列は対称になる）。\\( p(I) \\) は元素の周期で、H は 1、Z ≤ 10 は 2、Z ≤ 18 は 3、Z ≤ 36 は 4、Z ≤ 54 は 5、Z ≤ 86 は 6、Z ≤ 103 は 7、
  希ガス（He、Ne、Ar、Kr、Xe、Rn）は 0 とします。\\( B \\)（bohr）は対称な表です。

| \\( p \\) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| 1 | | −0.2573 | 0.3401 | 0.6937 | 0.7126 | 0.8335 | 0.9491 | 1.0 |
| 2 | | | 0.9652 | 1.2843 | 1.4625 | 1.6549 | 1.7190 | 2.0 |
| 3 | | | | 1.6925 | 1.8238 | 2.1164 | 2.3185 | 2.5 |
| 4 | | | | | 2.0203 | 2.2137 | 2.5206 | 2.7 |
| 5 | | | | | | 2.3718 | 2.5110 | 2.7 |
| 6 | | | | | | | 2.5200 | 2.7 |
| 7 | | | | | | | | 2.9 |

- 原子対の寄与は一様な並進 \\( \mathbf t \\) について \\( H\mathbf t = 0 \\) を満たし、0.1 の対角項だけが並進方向の曲率を与えます（`tools/check_geometry_md.py`）。
  同じ原子の周期像との対（\\( J = I \\)、\\( \mathbf L \ne 0 \\)）の寄与は対角ブロックで打ち消して 0 になります。
- 固定された成分（4.10 節）の行と列は 0 にし、対角を 1 にします。
- 原子対の範囲が PAO のカットオフ半径の和で決まるため、遠い対（\\( r \\) が 10 bohr 以上）も小さな寄与を持ちます。
  \\( r \le B \\) の対では \\( k \\) が負または発散します（通常の結合距離では起きない）。再実装では \\( r - B \\) に下限を設けます〔設計〕。
- `FF` という第 3 の選択肢（Lennard-Jones 型の力場を力に当てはめて Hessian を作る）がありますが、公開資料に記載がなく、この章では扱いません（〔未確定〕）。

### 4.10 拘束 `MD.Fixed.XYZ`〔公開〕〔規約〕

〔公開〕マニュアル 18.4 節：ブロックの各行は「原子番号、x、y、z のフラグ」で、1 が固定、0 が自由。既定はすべて自由です。
すべての構造最適化と MD で使えます。

〔規約〕

- 原子数ぶんの行を読みます（1 列目の番号は読み飛ばし、行の順で原子に割り当てる）。行が足りなければ入力エラーです。
- 構造最適化：固定成分の勾配を 0 として扱い、座標を動かしません。最大の力・最大変位の判定から除きます。
- MD：固定成分は位置を更新せず、速度を 0 に保ちます（ただし 6.1 節の例外）。
- `RFC5`–`RFC7` では、フラグは **分数座標の成分**（\\( \mathbf a\_1,\mathbf a\_2,\mathbf a\_3 \\) 方向）に対する固定として働きます（5.3 節）。

### 4.11 再開〔公開〕〔規約〕

〔公開〕マニュアル 17.3・18.5 節：各ステップで `System.Name.dat#`（最後の構造を入れた入力ファイル）が書かれ、
近似 Hessian も毎ステップ保存されて、`System.Name.dat#` から再開すると自動的に引き継がれます。分割して実行しても反復回数は変わりません。

〔規約〕構造最適化（`Opt`、`DIIS`、`BFGS`、`RF`、`EF`）の `.dat#` には `geoopt.restart on` が加わります。
`geoopt.restart on` のとき、最初のステップで `System.Name_rst/` に保存した最適化の状態（段階の別、段階内のステップ数、前因子、エネルギーの履歴、座標・勾配の履歴、Hessian）を読み込みます。
保存形式はこのノートで自由に定めます〔設計〕。

## 5. 格子の最適化

### 5.1 変数と勾配〔公開〕〔規約〕

〔公開〕マニュアル 19 節：格子の最適化は共線スピン（＋U を含む）でだけ使え、DFT-D2・D3 とは併用できません。
収束の判定は、原子座標または格子ベクトルに関するエネルギーの勾配の最大値と `MD.Opt.criterion` の比較です。
マニュアル 19.4 節：`MD.applied.pressure`（GPa、既定 0）を与えると、エンタルピー \\( H = U + P\Omega \\) を最小にします。

〔規約〕

- 格子ベクトルに関する勾配は [全エネルギーと力](energy-force.md) の式 \\( \partial U/\partial A\_{ij} = \sum\_k (A^{-1})\_{ki}\Sigma\_{kj} \\) です。
- 格子の最適化では応力の計算が自動的に有効になり、加えた圧力が応力に入ります：
  \\( \Sigma\_{\gamma\gamma} \mathrel{+}= P\Omega f\_\gamma \\)（\\( f\_\gamma \\) は `MD.applied.pressure.flag` の \\( \gamma \\) 番目、既定 `1 1 1`）。
  エネルギーには \\( U\_{pV} = P\Omega \\)（フラグによらない）を加え、`Enpy` = `Utot` + `UpV` を出力します（[全エネルギーと力](energy-force.md) の `UpV` の〔未確定〕はこれで解消：格子の最適化でだけ応力に入る）。
- 〔公開〕マニュアル 19.4 節：`MD.applied.pressure.flag` は、直方晶なら軸ごとに圧力をかけ、それ以外では他の 2 本の格子ベクトルに垂直な方向にかける、と説明されています。
  既存実装は上のとおりデカルトの対角成分 \\( \Sigma\_{\gamma\gamma} \\) に足すので、この説明が正確なのは格子ベクトルが座標軸に沿う場合です。
- 格子を動かすとき、分数座標 \\( \mathbf s\_I \\)（原点 \\( \mathbf O \\) からの）を保って \\( \mathbf R\_I = \mathbf O + A^{\mathsf T}\mathbf s\_I \\) で原子を動かします。
- 格子の最適化の各ステップでは、原子対の集合（近接リスト）を作り直します（通常の構造最適化・MD と同じ扱い）。

### 5.2 `OptC1`–`OptC7`（最急降下）〔公開〕〔規約〕

〔公開〕マニュアル 19.1 節の表：`OptC1` は格子に拘束なし・内部座標固定、`OptC2` は角度も固定し長さだけ、
`OptC3` はさらに \\( \lvert\mathbf a\_1\rvert = \lvert\mathbf a\_2\rvert = \lvert\mathbf a\_3\rvert \\)、`OptC4` は \\( \lvert\mathbf a\_1\rvert = \lvert\mathbf a\_2\rvert \ne \lvert\mathbf a\_3\rvert \\)、
`OptC5` は格子と内部座標を同時に、`OptC6` は \\( \mathbf a\_3 \\) 固定、`OptC7` は \\( \mathbf a\_2, \mathbf a\_3 \\) 固定。

〔規約〕\\( D\_{ij} = \partial H/\partial A\_{ij} \\)（固定成分は 0）、\\( s \\) を前因子として：

| 方式 | 収束判定の量 \\( G\_{\max} \\) | 更新 |
|---|---|---|
| `OptC1` | \\( \max\_{ij}\lvert D\_{ij}\rvert \\) | \\( A\_{ij} \leftarrow A\_{ij} - c s D\_{ij} \\)（自由な成分） |
| `OptC2` | \\( \max\_i \lvert d\_i\rvert \\)、\\( d\_i = (\mathbf a\_i\cdot\mathbf D\_i)/\lVert\mathbf a\_i\rVert \\) | \\( \mathbf a\_i \leftarrow (1 - c s d\_i)\mathbf a\_i \\)（3 成分とも自由な \\( \mathbf a\_i \\) だけ） |
| `OptC3` | 同上（固定フラグを見ない） | 同上の後、3 本の長さをその平均にそろえる |
| `OptC4` | 同上 | 同上の後、\\( \mathbf a\_1, \mathbf a\_2 \\) の長さをその平均にそろえる |
| `OptC5`–`OptC7` | \\( \max(\max\_{ij}\lvert D\_{ij}\rvert,\ \max\_{I\alpha}\lvert g\_{I\alpha}\rvert) \\) | 格子は `OptC1` と同じ。原子は新しい格子で \\( \mathbf s\_I \\) を保って動かした後、\\( \mathbf R\_I \leftarrow \mathbf R\_I - 0.2\\, c s\\, \mathbf g\_I \\) |

- `OptC6`、`OptC7` は `OptC5` に \\( \mathbf a\_3 \\)（と \\( \mathbf a\_2 \\)）の全成分の固定を加えたものです。
- 前因子 \\( s \\) は 4.2 節と同じ初期値・範囲ですが、エネルギーには \\( H = U + P\Omega \\) を使い、調整の規則が違います：
  \\( H\_n > H\_{n-1} \\) かつ \\( n \bmod 3 = 1 \\) なら \\( s \leftarrow s/5 \\)、\\( H\_{n-1} < H\_{n-2} \\) かつ \\( H\_n < H\_{n-1} \\) で \\( n \bmod 4 = 1 \\) なら \\( s \leftarrow 2.5 s \\)。
- \\( c \\) は \\( s G\_{\max} \\) が 0.5 を超えたときの縮小率 \\( 0.5/(sG\_{\max}) \\)（4.1 節と同じ考え方。格子の成分も bohr で比べる）。
- `OptC2` の更新は \\( \mathbf a\_i \\) 方向の単位ベクトルではなく \\( \mathbf a\_i \\) 自身に \\( c s d\_i \\) を掛けた量だけ動かすので、
  長さの変化は \\( c s d\_i\lVert\mathbf a\_i\rVert \\) です（次元が整合しないが、既存実装の規約）。
- `OptC3`・`OptC4` で長さをそろえる操作は、初期構造が条件を満たしていなくても行われます（入力の検査はない）。

### 5.3 `RFC5`–`RFC7`〔公開〕〔規約〕

〔公開〕マニュアル 19.1 節：格子と内部座標を同時に、RF と DIIS の混成法（Hessian は BFGS で更新）で最適化します。
任意の拘束には `RFC5` と `MD.Fixed.Cell.Vectors` の組み合わせを勧めています。

〔規約〕

- 変数は分数座標 \\( \mathbf s\_I \\)（\\( 3N \\) 個）と、動かす格子ベクトル（`RFC5` は 3 本、`RFC6` は \\( \mathbf a\_1,\mathbf a\_2 \\)、`RFC7` は \\( \mathbf a\_1 \\)）のデカルト成分です。
  勾配は \\( \partial H/\partial s\_{Ii} = \sum\_\alpha g\_{I\alpha} A\_{i\alpha} \\) と \\( D\_{ij} \\) です。
- 更新は 4.3 節の GDIIS と 4.6 節の RF ステップを、この変数の組に対して行います（BFGS の分母の扱いも 4.6 節と同じ）。
  最大変位の判定は、更新後のデカルト座標と格子成分の変化の最大（bohr）で 0.5 と比べます。
- 収束判定の量は \\( \max(\max\_{ij}\lvert D\_{ij}\rvert,\ \max\_{I\alpha}\lvert g\_{I\alpha}\rvert) \\) です（デカルトの勾配で判定する）。
- 準 Newton 段階に入る前（4.8 節）は `OptC5` 型の最急降下を使います。
- 初期 Hessian（準 Newton 段階の最初のステップでだけ作る）：
  - `iden`：単位行列（\\( 3N + 9 \\) 次元）。
  - `Schlegel`：原子の部分は対角 0.1、格子の部分は対角 0.05 から始め、4.9 節の原子対の寄与を加えた後、原子の部分を
    分数座標に変換し（\\( H^{(s)}\_{Ii,Jj} = \sum\_{\alpha\beta} A\_{i\alpha}H\_{I\alpha,J\beta}A\_{j\beta} \\)）、
    格子どうし・原子と格子の間の部分を、同じ結合伸縮のモデルをひずみで 2 回微分した量から作ります。
    この交差項の詳細は収束の経路にしか影響しないので、この章では定めません（〔未確定〕）。
- 固定された分数座標成分は、Hessian の行と列を 0 にし、対角を 200 にします。
- `MD.Fixed.Cell.Vectors` で固定した格子成分は勾配を 0 にします。

### 5.4 `MD.Fixed.Cell.Vectors`〔公開〕〔規約〕

〔公開〕マニュアル 19.3 節：3 行 3 列のフラグで、行 \\( i \\) が \\( \mathbf a\_i \\)、列が x、y、z 成分。1 が固定、0 が自由。
`OptC1`、`OptC2`、`OptC5`、`RFC5` で使えます。

〔規約〕このブロックがあれば、`OptC6`・`OptC7`・`RFC6`・`RFC7` が決める固定より優先します。
`OptC2`–`OptC4` では、ベクトル単位の更新なので、3 成分とも自由なベクトルだけが動きます（`OptC3`・`OptC4` の長さをそろえる操作はフラグを見ない）。

### 5.5 `EvsLC` と `DF`〔規約〕

- `EvsLC`：各ステップの後、最初の格子 \\( A^{(0)} \\) を使って \\( \mathbf a\_i \leftarrow \mathbf a\_i \pm (\delta/100)\\,\mathbf a\_i^{(0)} \\) とします。
  \\( \delta \\) は `MD.EvsLC.Step`（%、既定 0.4）、符号は `MD.EvsLC.flag`（3 個の整数、1 で伸ばす、−1 で縮める、0 で固定。既定 `1 1 1`）。
  分数座標は保ちます。各ステップの格子ベクトル（Å）と `Utot` を `System.Name.EvsLC` に追記します。
- `DF`（デルタ因子用）：`MD.maxIter` は 7 でなければ入力エラーです。ステップ 1 は入力の格子、ステップ 2–7 は
  体積を 0.94、0.96、0.98、1.02、1.04、1.06 倍（格子を体積比の 1/3 乗で等方的に拡大縮小）した格子で計算し、
  各ステップの原子あたり体積（Å³）と原子あたり `Utot`（eV）を `System.Name.DF` に追記します。

## 6. 分子動力学

### 6.1 初速度〔公開〕〔規約〕

〔公開〕マニュアル 20.8 節：`MD.Init.Velocity` ブロックの各行は「原子番号、\\( v\_x, v\_y, v\_z \\)（m/s）」で、`MD.Fixed.XYZ` と併用できます。

〔規約〕

- `MD.Init.Velocity` があれば、それを現在の速度とします（固定成分は 0 にする）。
  `MD.Init.Velocity.Prev`（同じ形式）は 1 ステップ前（半ステップ）の速度で、再開ファイルが書き込みます。
- どちらもなく、`MD.Type` が熱浴付き（`NVT_*`、`NPT_*`）なら、乱数で初速度を作ります：
  1. 各原子に \\( \lVert\mathbf v\_I\rVert = \sqrt{3k\_BT\_1/m\_I} \\)（\\( T\_1 \\) は `MD.TempControl` の最初の行の温度）、向きは \\( [-\tfrac12,\tfrac12]^3 \\) の一様乱数ベクトルを規格化したもの。
  2. 全運動量を 0 にする：\\( \mathbf v\_I \leftarrow \mathbf v\_I - \bar{\mathbf P}/m\_I \\)、\\( \bar{\mathbf P} = \frac1N\sum\_J m\_J\mathbf v\_J \\)。
  3. \\( T\_A \\)（2 節。固定成分を除いた \\( K \\)）がちょうど \\( T\_1 \\) になるよう全体を定数倍する。
  - 乱数列は C ライブラリの擬似乱数（種 1）で、環境に依存します。照合には `MD.Init.Velocity` を明示した入力を使います〔設計〕。
  - この場合、固定成分の速度は 0 にされません（\\( K \\) には入らないが、`NVT_VS` 系では温度に入る。11 節）。
- それ以外（`NVE` など）では初速度 0 です。

### 6.2 `NVE`（速度 Verlet）〔公開〕〔規約〕

〔公開〕マニュアル 20.1 節：微視的正準（NVE）の MD。`.ene` の 14–17 列（`Ukc`、`Utot`、`Utot+Ukc`、Fermi 準位）で保存を確かめます。

〔規約〕速度 Verlet 法です。\\( \mathbf a\_n = -\mathbf g\_n/m \\)（1 節の補正後）として、ステップ \\( n \\) で

$$
\mathbf v\_n = \mathbf v\_{n-1/2} + \tfrac12\mathbf a\_n\Delta t\ \ (n\ge 2),\qquad
\mathbf v\_{n+1/2} = \mathbf v\_n + \tfrac12\mathbf a\_n\Delta t,\qquad
\mathbf x\_{n+1} = \mathbf x\_n + \mathbf v\_{n+1/2}\Delta t
$$

- \\( n = 1 \\) では入力の速度を \\( \mathbf v\_1 \\) とします。\\( K \\) と温度 \\( T\_A \\) は \\( \mathbf v\_n \\) で求め、`.ene` に書きます。
- 保持する速度は、ステップの終わりで \\( \mathbf v\_{n+1/2} \\) です（次のステップの `.md` にはこの半ステップの速度が出る）。
- 固定成分は位置・速度とも更新しません。
- 調和振動子でエネルギー誤差が \\( O(\Delta t^2) \\) で有界であることを `tools/check_geometry_md.py` で確認しています。

### 6.3 温度制御表 `MD.TempControl`〔公開〕〔規約〕

〔公開〕マニュアル 20.2・20.3 節。ブロックの最初の数は行数 \\( L \\) で、続く行の形式が方式で異なります。

- 速度スケーリング（`NVT_VS`、`NVT_VS2`、`NVT_VS4`、`NPT_VS_PR`、`NPT_VS_WV`）：各行「\\( N\_i \\)（区間の終わりのステップ）、\\( I\_i \\)（スケーリングの間隔）、\\( T\_i \\)（K）、\\( \alpha\_i \\)」。
- Nosé–Hoover・Langevin（`NVT_NH`、`NVT_Langevin`、`NPT_NH_PR`、`NPT_NH_WV`）：各行「\\( S\_i \\)（ステップ）、\\( T\_i \\)（K）」。
- いずれも行の間の目標温度はステップについて線形補間します。

〔規約〕補間の区間の取り方が方式で違います。\\( (N\_0, T\_0) = (0, T\_1) \\)（または \\( (S\_0, T\_0) = (0, T\_1) \\)）を補います。

| 方式 | ステップ \\( n \\) が属する区間 | 区間の外 |
|---|---|---|
| 速度スケーリング（NVT・NPT） | \\( N\_{i-1} < n \le N\_i \\) | スケーリングしない（`NVT_VS` 系は 6.4 節の積分だけ、マニュアルの言う NVE への移行） |
| `NVT_NH` | \\( S\_i \le n < S\_{i+1} \\)（\\( i = 1,\dots,L-1 \\)） | 熱浴を切って速度 Verlet（6.2 節）で進める |
| `NVT_Langevin` | 同上 | 目標温度は直前に設定された値のまま、熱浴は働き続ける |
| `NPT_NH_*` | \\( S\_{i-1} < n \le S\_i \\) | 目標温度は直前の値のまま、熱浴は働き続ける |

### 6.4 速度スケーリング `NVT_VS`、`NVT_VS2`、`NVT_VS4`〔公開〕〔規約〕

〔公開〕マニュアル 20.2 節：区間 \\( i \\) では \\( I\_i \\) ステップごとに速度を \\( s \\) 倍し、
スケーリング後の温度を \\( T\_0 + (T - T\_0)\alpha\_i \\) にします（\\( \alpha = 0 \\) で目標温度に合わせ、\\( \alpha = 1 \\) で何もしない）：

$$
s = \sqrt{\frac{T\_0 + (T - T\_0)\alpha\_i}{T}}
$$

（マニュアルの HTML の式は平方根が崩れて読めますが、温度をこの値にする係数は平方根付きです。`tools/check_geometry_md.py`。）
`NVT_VS2` は原子ごとの温度で原子ごとに、`NVT_VS4`（マニュアル 20.4 節）は原子群ごとの温度で群ごとにスケーリングします。
原子群は `MD.num.AtomGroup`（群の数）と `MD.AtomGroup`（各行「原子番号、群番号」、群番号は 1 から連番）で指定し、
すべての群が同じ `MD.TempControl` に従います。

〔規約〕

- 積分の形（3 方式共通）：\\( \mathbf u \\) を保持する速度として
  - \\( n = 1 \\)：\\( \mathbf x\_2 = \mathbf x\_1 + \mathbf u\_1\Delta t + \tfrac12\mathbf a\_1\Delta t^2 \\)、\\( \mathbf u \leftarrow \mathbf u\_1 + \mathbf a\_1\Delta t \\)。スケーリングなし。
  - \\( n \ge 2 \\)：\\( \mathbf u \leftarrow \mathbf u + \mathbf a\_n\Delta t \\) → 速度の平均を引く（下記）→ 温度を計算 → スケーリングの判定と \\( \mathbf u \leftarrow s\mathbf u \\) →
    \\( \mathbf x\_{n+1} = \mathbf x\_n + \mathbf u\Delta t + \tfrac12\mathbf a\_n\Delta t^2 \\)。
  - これは速度 Verlet では **ありません**。スケーリングがない場合、\\( \mathbf x\_{n+1} - 2\mathbf x\_n + \mathbf x\_{n-1} = \tfrac12(3\mathbf a\_n - \mathbf a\_{n-1})\Delta t^2 \\) となり、
    調和振動子 \\( \omega \\) では特性根の絶対値が \\( \sqrt{1 - (\omega\Delta t)^2/2} \\)、つまり 1 ステップごとにエネルギーが \\( (\omega\Delta t)^2/2 \\) の割合で失われます
    （周期 33 fs の振動、\\( \Delta t = 0.5 \\) fs で 1 ステップ 0.45%）。速度スケーリングがこの散逸を補っています。`tools/check_geometry_md.py` で確認しています。
- 速度の平均：自由な成分について、デカルト成分ごとに速度の単純平均（質量の重みなし）を引きます。全運動量 0 とは一致しません（質量が異なる場合）。
- 温度：`NVT_VS` は \\( T\_A \\)、ただし \\( K \\) は **全成分**（固定成分を含む）で和を取ります。
  `NVT_VS2` は原子ごとの \\( T\_I = K\_I/(\tfrac32k\_B) \\)（平均を引く前の速度で計算）、`NVT_VS4` は群ごとの \\( T\_g = K\_g/(\tfrac32 k\_B N\_g) \\)（平均を引いた後）。
- スケーリングは、区間 \\( i \\) で \\( (n - N\_{i-1}) \bmod I\_i = 0 \\) のステップだけ行います。目標温度はそのステップの補間値です。
- 最終ステップ（\\( n = \\) `MD.maxIter`）では、スケーリングと位置の更新をしません。

### 6.5 Nosé–Hoover `NVT_NH`〔公開〕〔規約〕

〔公開〕マニュアル 20.3 節：熱浴の質量 `NH.Mass.HeatBath`（既定 20.0、単位は 統一原子質量単位 × bohr²）。拡張系の Hamiltonian が保存量として `.ene` に出ます。
運動方程式は Nosé（Mol. Phys. **52**, 255 (1984)）、Hoover（Phys. Rev. A **31**, 1695 (1985)）の形

$$
\dot{\mathbf R}\_I = \mathbf v\_I,\qquad
\dot{\mathbf v}\_I = \frac{\mathbf F\_I}{m\_I} - \zeta\mathbf v\_I,\qquad
\dot\zeta = \frac{2K - gk\_BT\_0}{Q},\qquad
\dot\eta = \zeta
$$

で、\\( H' = U + K + \tfrac12 Q\zeta^2 + gk\_BT\_0\eta \\) が保存されます（\\( T\_0 \\) 一定のとき）。

〔規約〕

- \\( g = 3N \\)（固定成分があっても全原子数で数える）、\\( Q \\) は `NH.Mass.HeatBath` を \\( m\_e \\)·bohr² に換算した値です。
- 離散化（速度 Verlet 型）。ステップ \\( n \\) の力が得られた時点で：
  1. **補正**（\\( n \ge 2 \\)）：前のステップの予測値 \\( \zeta' \\) と半ステップ速度 \\( \mathbf v\_{n-1/2} \\) から、連立方程式
     \\( \mathbf v\_n = \mathbf v\_{n-1/2} + \tfrac12\Delta t\\,(\mathbf a\_n - \zeta\_n\mathbf v\_n) \\)、
     \\( \zeta\_n = \zeta' + \Delta t\\,(K(\mathbf v\_n) - \tfrac12 gk\_BT\_0)/Q \\)
     を Newton 法で解きます（反復は最大 21 回。収束の判定の細部はビット一致に関わるだけなので定めない）。
  2. 速度の平均を引き（6.4 節と同じ）、\\( K \\)、\\( T\_A \\)、\\( H' \\) を計算して `.ene` に書く。
  3. **予測**（最終ステップ以外）：\\( \mathbf v\_{n+1/2} = \mathbf v\_n + \tfrac12\Delta t\\,(\mathbf a\_n - \zeta\_n\mathbf v\_n) \\)、\\( \mathbf x\_{n+1} = \mathbf x\_n + \mathbf v\_{n+1/2}\Delta t \\)、
     \\( \zeta' = \zeta\_n + \Delta t\\,(K(\mathbf v\_n) - \tfrac12 gk\_BT\_0)/Q \\)、\\( \eta\_{n+1} = \eta\_n + \zeta'\Delta t \\)。
  - 1 と 3 を合わせると \\( \zeta\_{n+1} = \zeta\_n + \tfrac{\Delta t}{2}[\dot\zeta\_n + \dot\zeta\_{n+1}] \\)（台形則）で、上の連続方程式の \\( \dot\zeta = (2K - gk\_BT\_0)/Q \\) に一致します。
    `tools/check_geometry_md.py` で、この離散化が \\( H' \\) を \\( O(\Delta t^2) \\) で保存することを確認しています。
- \\( n = 1 \\) では \\( \zeta = \eta = 0 \\) から始めます（入力の `NH.czeta` 等は再開のとき、つまり通算 \\( n \ge 2 \\) でだけ使われる）。
- 区間の外（6.3 節）では \\( \zeta \\) を使わず速度 Verlet で進めます。`.ene` の \\( H' \\) はそのときの \\( \zeta, \eta \\) で計算されます。
- `.ene` の 20・21・22 列は \\( \zeta \\)、予測値 \\( \zeta' \\)、\\( H' \\) です。

### 6.6 Langevin `NVT_Langevin`〔規約〕

公開資料に説明がありません（〔規約〕のみ）。運動方程式は標準的な Langevin 方程式

$$
m\_I\dot{\mathbf v}\_I = \mathbf F\_I - m\_I\gamma\mathbf v\_I + \mathbf R\_I(t),\qquad
\langle R\_{I\alpha}(t)R\_{J\beta}(t')\rangle = 2m\_I\gamma k\_BT\_0\\,\delta\_{IJ}\delta\_{\alpha\beta}\delta(t-t')
$$

です。\\( \gamma \\) は `Langevin.Friction.Factor`（原子単位の時間の逆数、既定 0.001、すなわち約 41 ps⁻¹）。

- 各ステップで、自由な成分に標準偏差 \\( \sqrt{2m\_I\gamma k\_BT\_0/\Delta t} \\) の Gauss 乱数（Box–Muller の片側）を力に加えます。以下 \\( \tilde{\mathbf a} \\) はこれを含む加速度です。
- 離散化：
  - \\( n = 1 \\)：\\( \mathbf x\_2 = \mathbf x\_1 + \mathbf u\Delta t + \tfrac12\Delta t^2(\tilde{\mathbf a}\_1 - \gamma\mathbf u) \\)、\\( \mathbf u \leftarrow \mathbf u + \tilde{\mathbf a}\_1\Delta t \\)。
  - \\( n \ge 2 \\)：\\( \mathbf u \leftarrow \mathbf u + \tfrac12(\tilde{\mathbf a}\_n + \tilde{\mathbf a}\_{n-1})\Delta t - \gamma(\mathbf x\_n - \mathbf x\_{n-1}) \\)、
    \\( \mathbf x\_{n+1} = \mathbf x\_n + (\mathbf u\Delta t + \tfrac12\tilde{\mathbf a}\_n\Delta t^2)/(1 + \tfrac12\gamma\Delta t) \\)。
  - 最終ステップでも位置を更新します（他の方式と異なる）。速度の平均は引きません。
- \\( K \\)、\\( T\_A \\) は更新後の \\( \mathbf u \\) で計算します。
- 乱数は 6.1 節と同じ C ライブラリの擬似乱数列の続きです。軌跡の照合はできないので、統計量（平均温度など）で比べます〔設計〕。

### 6.7 NPT `NPT_VS_PR`、`NPT_VS_WV`、`NPT_NH_PR`、`NPT_NH_WV`〔公開〕〔規約〕

〔公開〕マニュアル 20.6・20.7 節：熱浴は速度スケーリング（VS）または Nosé–Hoover（NH）、
圧力浴は Parrinello–Rahman（PR、J. Appl. Phys. **52**, 7182 (1981)）または Wentzcovitch（WV、Phys. Rev. B **44**, 2358 (1991)）。
キーワード：

| キーワード | 既定 | 意味 |
|---|---|---|
| `MD.TempControl` | — | 6.3 節 |
| `MD.CellPressureControl` | — | 行数の後に「ステップ、圧力（GPa）」の行。間は線形補間。正の圧力で圧縮 |
| `NPT.Mass.Barostat` | 自動 | バロスタットの質量 \\( W \\)（統一原子質量単位） |
| `NPT.LatticeRestriction` | `none` | `none`、`cubic`、`ortho` |
| `MD.TempTolerance` | 100 K | \\( \lvert T - T\_0\rvert \\) がこれ未満ならスケーリングしない（VS） |
| `MD.applied.pressure.flag` | `1 1 1` | 圧力をかける成分 |
| `NH.Mass.HeatBath` | 20 | NH の熱浴の質量（6.5 節と同じ単位） |

拘束（`MD.Fixed.XYZ`）は NPT と併用できません（マニュアル 20.6 節）。

〔規約〕変数は分数座標 \\( \mathbf s\_I \\)、その時間微分 \\( \dot{\mathbf s}\_I = A^{-\mathsf T}\mathbf v\_I \\)、格子 \\( A \\) とその速度 \\( \dot A \\) です
（格子の速度は再開ファイルの `Atoms.Unitvectors.Velocity` ブロックに保存され、最初は 0）。
\\( G = AA^{\mathsf T} \\)（\\( G\_{ij} = \mathbf a\_i\cdot\mathbf a\_j \\)）として、原子の運動方程式は

$$
\ddot{\mathbf s}\_I = -\frac{A^{-\mathsf T}\mathbf g\_I}{m\_I} - G^{-1}\dot G\\,\dot{\mathbf s}\_I\ \ [-\ \zeta\dot{\mathbf s}\_I\ \text{(NH)}]
$$

格子の運動方程式（PR）は、内部圧力テンソル

$$
\Pi\_{\alpha\beta} = \frac1\Omega\sum\_I m\_I v\_{I\alpha}v\_{I\beta} - \frac{\Sigma\_{\alpha\beta}}{\Omega},\qquad \mathbf v\_I = A^{\mathsf T}\dot{\mathbf s}\_I
$$

を使って

$$
W\ddot A\_{i\alpha} = \sum\_\beta\left(\Pi\_{\alpha\beta} - P f\_\alpha\delta\_{\alpha\beta}\right)\sigma\_{i\beta}
$$

です（\\( \sigma\_{i\beta} \\) は \\( S \\) の要素。\\( \partial\Omega/\partial A\_{i\beta} = \sigma\_{i\beta} \\)、`tools/check_geometry_md.py`）。
`cubic` では右辺の括弧の行列の非対角を 0 にし対角をその平均に、`ortho` では非対角を 0 にします。

- 時間積分（1 ステップ。半陰的 Euler）：
  1. \\( \ddot{\mathbf s} \\) を現在の量で計算し、\\( \dot{\mathbf s} \leftarrow \dot{\mathbf s} + \ddot{\mathbf s}\Delta t \\)。
  2. \\( K \\) と \\( T\_A \\) を新しい \\( \dot{\mathbf s} \\) と現在の \\( A \\) で計算。
     VS：区間内で \\( n = 1 \\) か、\\( (n - N\_{i-1}) \bmod I\_i = 0 \\) かつ \\( \lvert T\_0 - T\rvert > \\) `MD.TempTolerance` のとき \\( \dot{\mathbf s} \\) を \\( s \\) 倍（\\( T < 10^{-8}T\_0 \\) なら倍率 1）。
     NH：\\( \zeta \leftarrow \zeta + Nk\_B(T\_A - T\_0)\Delta t/Q \\)、\\( \eta \leftarrow \eta + \zeta\Delta t \\)（\\( n = 1 \\) で \\( \zeta = \eta = 0 \\) から）。
  3. \\( \mathbf s \leftarrow \mathbf s + \dot{\mathbf s}\Delta t \\)。
  4. \\( \ddot A \\) を新しい \\( \dot{\mathbf s} \\)、現在の \\( A \\) と応力で計算し、\\( \dot A \leftarrow \dot A + \ddot A\Delta t \\)、\\( A \leftarrow A + \dot A\Delta t \\)。
  5. \\( \mathbf R\_I = \mathbf O + A^{\mathsf T}\mathbf s\_I \\)、\\( \mathbf v\_I = A^{\mathsf T}\dot{\mathbf s}\_I \\)。
- 外圧 \\( P \\)：`MD.CellPressureControl` の行 \\( (S\_i, P\_i) \\) に \\( (0, P\_1) \\) を補い、\\( S\_{i-1} < n \le S\_i \\) で補間します。
  最後の行のステップを超えた場合の値は定義されていません（〔規約〕、再実装では最後の値を保つ〔設計〕）。`MD.applied.pressure` は NPT では使いません。
- NH の熱浴は \\( \dot\zeta = Nk\_B(T\_A - T\_0)/Q = (2K - 3Nk\_BT\_0)/(3Q) \\) です。これは 6.5 節の形で質量を \\( 3Q \\) にしたものに当たり、
  保存量は \\( \tfrac12(3Q)\zeta^2 \\) を含みます。既存実装が出力する値は \\( \tfrac12 Q\zeta^2 + 3Nk\_BT\_0\eta \\) を使うので、厳密には保存しません（`tools/check_geometry_md.py` の 5 番）。
- 既定の \\( W \\)：PR では \\( W = \frac{3}{4\pi^2}\sum\_I m\_I \\)、WV では \\( W = \frac{3}{4\pi^2}\sum\_I m\_I\times\frac{0.01}{\Omega^{2/3}} \\)（\\( \Omega \\) は最初の体積、bohr³）。
  〔公開〕マニュアル 20.6 節の Eq. 27 は PR の式と一致します。20.7 節の WV の既定値の式は表示が崩れていて読み取れません。
- WV の格子の運動方程式は、格子の運動エネルギーを \\( \tfrac12 W\sum\_{ij}(\dot{\mathbf a}\_i\cdot\dot{\mathbf a}\_j)F\_{ij} \\)、\\( F = SS^{\mathsf T} \\)（現在の格子）とした Lagrangian の Euler–Lagrange 方程式

$$
W\sum\_j\left(F\_{ij}\ddot{\mathbf a}\_j + \dot F\_{ij}\dot{\mathbf a}\_j\right) = \mathbf f\_i + \tfrac12 W\sum\_{jk}(\dot{\mathbf a}\_j\cdot\dot{\mathbf a}\_k)\\,\frac{\partial F\_{jk}}{\partial\mathbf a\_i},\qquad
f\_{i\alpha} = \sum\_\beta\left(\Pi\_{\alpha\beta} - Pf\_\alpha\delta\_{\alpha\beta}\right)\sigma\_{i\beta}
$$

  を、\\( F^{-1} = G/\Omega^2 \\) を掛けて \\( \ddot{\mathbf a}\_i \\) について解いたものです（\\( \mathbf f\_i \\) は PR の右辺と同じ。
  \\( \partial\boldsymbol\sigma\_j/\partial a\_{i\alpha} \\) は \\( (\sigma\_{i\alpha}\boldsymbol\sigma\_j - \sigma\_{j\alpha}\boldsymbol\sigma\_i)/\Omega \\)）。
  格子の制限（`cubic`、`ortho`）は \\( \mathbf f \\) の括弧の行列にかけます。〔規約〕既存実装では \\( \dot F\_{ij} \\) を
  \\( 2\dot{\boldsymbol\sigma}\_i\cdot\boldsymbol\sigma\_j \\)（対称化していない）で計算しています。
- 出力：`System.Name.npt.ene` に 1 行ずつ「ステップ、時刻 \\( n\Delta t \\)（fs）、体積（bohr³）、\\( T\_A \\)（K）、内部圧力 \\( -\mathrm{tr}\\,\Sigma/(3\Omega) \\)（GPa、係数 29.42×10³）、
  `Utot`、`Ukc`、`Ubaro`、`Wbaro` = \\( P\Omega \\)、その和（NH では熱浴のエネルギーも加える）」（Hartree）。
  PR の `Ubaro` は \\( \tfrac12 W\sum\_{i}\lVert\dot{\mathbf a}\_i\rVert^2 \\)。WV の `Ubaro` は最初の格子の \\( F^{(0)} \\) から作った量で、上の運動エネルギーとは一致しません（〔未確定〕）。
  見出しの 2 行（方式名、時間刻み、\\( W \\)（統一原子質量単位））がファイルの先頭に付きます。
- NPT では `.ene` を書きません。`.md` の速度の列は 0 です（10 節）。

### 6.8 MD の再開〔規約〕

`System.Name.dat#` に次が書かれ、そのまま入力にすると続きから計算します。

- `MD.Current.Iter`（そのステップ番号）、`<Atoms.SpeciesAndCoordinates`（現在の座標、入力と同じ単位）、`scf.restart on`、`scf.fixed.grid`（原点）。
- `<Atoms.UnitVectors`（現在の格子）、NPT では `<Atoms.Unitvectors.Velocity`（\\( \dot A \\)）と、WV では `<NPT.WV.F0`。
- MD では `<MD.Init.Velocity`（保持している速度、m/s）。`NVE`、`NVT_NH`、`NVT_Langevin`、`NPT_NH_*` ではさらに `<MD.Init.Velocity.Prev`（半ステップの速度）。
  NH では `NH.R`（\\( \eta \\)）、`NH.nzeta`（\\( \zeta' \\)）、`NH.czeta`（\\( \zeta \\)）。

## 7. 周期境界での原子の折り返し〔規約〕

- 分子動力学（`NVE`、`NVT_*`、`NPT_*`）の各ステップの更新後、原点 \\( \mathbf O \\) からの分数座標
  \\( s\_i = (\mathbf R - \mathbf O)\cdot\mathbf b\_i/(2\pi) \\)（\\( \mathbf b\_i \\) は逆格子ベクトル）が範囲外の原子を戻します：
  \\( s\_i > 1 \\) なら \\( s\_i \leftarrow s\_i - \lfloor s\_i\rfloor \\)、\\( s\_i < 0 \\) なら \\( s\_i \leftarrow s\_i + \lfloor -s\_i\rfloor + 1 \\)。
  \\( s\_i = 1 \\) ちょうどはそのまま残ります（\\( s\_i = -1 \\) ちょうどは 1 になる）。
- 戻した原子については、密度の外挿に使う座標の履歴と、GDIIS の座標の履歴も同じ格子ベクトルだけずらし、外挿や履歴が連続になるようにします
  （[密度行列と電荷ミキシング](density-mixing.md) 7 節）。戻したかどうかは次のステップの SCF に伝えられ、前のステップの量からの再開のしかたに影響します（詳細は〔未確定〕）。
- 構造最適化と格子の最適化では折り返しません。NEB でも折り返しません（8 節）。
- 孤立系（`Cluster` など）でも同じ処理が箱に対して行われます。分子が箱の境界をまたぐと分子が分断されるので、再実装では孤立系では折り返さないことを推奨します〔設計〕。
- 入力で分数座標（`FRAC`）を使った場合、読み込み時にも \\( [0, 1) \\) に戻します（同じ規則）。

## 8. NEB

### 8.1 入力〔公開〕〔規約〕

〔公開〕マニュアル 51 節：反応物と生成物の構造をそれぞれ最適化してから、`MD.Type NEB` で両者をつなぐ最小エネルギー経路を求めます。
同じ単位胞・カットオフ・基底・擬ポテンシャル・電子温度を使います。周期系では座標を `Ang` か `AU` で与え、`FRAC` を使わないこと（像が不連続になる）。

| キーワード | 既定（既存実装） | マニュアルの記載 | 意味 |
|---|---|---|---|
| `MD.NEB.Number.Images` | 10 | 10 | 端を除いた像の数 \\( P \\) |
| `MD.NEB.Spring.Const` | 0.4 | 0.1 | ばね定数 \\( k \\)（Hartree/bohr²） |
| `MD.NEB.Parallel.Number` | 0（自動） | — | 像の並列化の群の数（8.4 節） |
| `MD.Opt.DIIS.History` | 2〔未確定〕 | 7 | DIIS の履歴 |
| `MD.Opt.StartDIIS` | 1 | 5 | このステップから DIIS＋BFGS |
| `MD.Opt.criterion` | 0.0003 | 0.0003 | Hartree/bohr |
| `MD.maxIter` | 1 | 1 | 最適化の反復の上限 |
| `MD.Fixed.XYZ` | なし | — | 固定する成分 |

〔規約〕

- 反応物は `<Atoms.SpeciesAndCoordinates`、生成物は `<NEB.Atoms.SpeciesAndCoordinates`（同じ形式）で与えます。
- 初期経路は、デカルト座標での直線補間 \\( \mathbf x\_p = \mathbf x\_0 + \frac{p}{P+1}(\mathbf x\_{P+1} - \mathbf x\_0) \\)（\\( p = 1,\dots,P \\)）です。
  `FRAC` の端点は読み込み時に \\( [0,1) \\) に戻されるので、周期の境界をまたぐ移動は不連続になります（マニュアルの注意の理由）。
  最小像規約は使いません。
- 〔公開〕マニュアル 51.5 節：`scf.restart on` のとき、`<NEB1.Atoms.SpeciesAndCoordinates` … `<NEBP.Atoms.SpeciesAndCoordinates` で各像の座標を与えられます。
- 各像（と、最初の 1 回だけ両端）は独立した通常の SCF 計算として実行され、`System.Name.dat_p`（入力）と `System.Name_p.out` が作られます（\\( p = 0,\dots,P+1 \\)）。

### 8.2 接線と力〔公開〕〔規約〕

〔公開〕NEB の技術ノート（T. Ozaki, 2011）とマニュアル 51.1 節：G. Henkelman and H. Jónsson, J. Chem. Phys. **113**, 9978 (2000) の
改良接線（Eq. 8–11）、垂直な力（Eq. 4）、ばねの平行な力（Eq. 12）を使います。climbing image はありません。

像 \\( p \\) のエネルギーを \\( V\_p \\)、\\( \boldsymbol\tau^+ = \mathbf x\_{p+1} - \mathbf x\_p \\)、\\( \boldsymbol\tau^- = \mathbf x\_p - \mathbf x\_{p-1} \\) として

$$
\boldsymbol\tau\_p = \begin{cases}
\boldsymbol\tau^+ & V\_{p-1} < V\_p < V\_{p+1} \\\\
\boldsymbol\tau^- & V\_{p+1} < V\_p < V\_{p-1} \\\\
\Delta V^{\max}\boldsymbol\tau^+ + \Delta V^{\min}\boldsymbol\tau^- & \text{極値で } V\_{p-1} \le V\_{p+1} \\\\
\Delta V^{\min}\boldsymbol\tau^+ + \Delta V^{\max}\boldsymbol\tau^- & \text{極値で } V\_{p-1} > V\_{p+1}
\end{cases}
$$

（\\( \Delta V^{\max/\min} = \max/\min(\lvert V\_{p+1} - V\_p\rvert, \lvert V\_{p-1} - V\_p\rvert) \\)）を規格化して \\( \hat{\boldsymbol\tau}\_p \\) とし、
最適化に使う勾配を

$$
\mathbf G\_p = \mathbf g\_p - (\mathbf g\_p\cdot\hat{\boldsymbol\tau}\_p)\hat{\boldsymbol\tau}\_p - k\left(\lVert\boldsymbol\tau^+\rVert - \lVert\boldsymbol\tau^-\rVert\right)\hat{\boldsymbol\tau}\_p
$$

とします（\\( \mathbf F^{\mathrm{NEB}}\_p = -\mathbf G\_p \\)）。

〔規約〕

- 「極値」は \\( V\_{p+1} \ge V\_p \le V\_{p-1} \\) または \\( V\_{p+1} < V\_p > V\_{p-1} \\) の場合です。
  \\( V\_p = V\_{p+1} > V\_{p-1} \\) のような等号の組み合わせはどの場合にも当たらず、接線が 0 になって規格化で発散します（11 節）。
- 固定成分の \\( \mathbf G\_p \\) は 0 にします。勾配は像ごとの SCF の勾配そのもので、総和 0 の補正はしません。

### 8.3 最適化〔公開〕〔規約〕

〔公開〕マニュアル 51.3 節、技術ノート：全像の座標をまとめた \\( 3NP \\) 次元のベクトルに対する DIIS と BFGS の混成法です。

〔規約〕

- 収束の判定：\\( \max\_{p,I}\lVert\mathbf G\_{p,I}\rVert < \\) `MD.Opt.criterion`（原子ごとのベクトルの大きさ）。
- \\( n < \\) `MD.Opt.StartDIIS`：最急降下 \\( \mathbf x \leftarrow \mathbf x - s\mathbf G \\)。前因子は基準 \\( s\_0 = (4\ \mathrm{fs})^2/1836.1526 \\) の \\( [0.005s\_0, 10s\_0] \\) に制限し、
  最初は \\( s = 1/(5\times0.529177249) \\)。\\( \lVert\mathbf G\rVert \\)（全体の Euclid ノルム）が前より増え \\( n \bmod 4 = 1 \\) なら \\( s/4 \\)、
  2 回続けて減り \\( n \bmod 4 = 1 \\) なら \\( 1.2s \\)。最大変位の制限はありません。
- \\( n \ge \\) `MD.Opt.StartDIIS`：4.3 節の GDIIS（\\( \mathbf g \\) の代わりに \\( \mathbf G \\)、全像をまとめる）で \\( \tilde{\mathbf x} \\)、\\( \tilde{\mathbf G} \\) を作り、
  4.5 節の BFGS の逆 Hessian（`MD.Opt.StartDIIS` のステップで \\( I \\)、曲率条件の確認なし）で
  \\( \mathbf x \leftarrow \tilde{\mathbf x} - f\\, H^{-1}\tilde{\mathbf G} \\)。\\( f \\) は最大の力 \\( F\_{\max} \\) で決めます：

| \\( F\_{\max} \\) | \\( >0.7 \\) | \\( >0.3 \\) | \\( >0.1 \\) | \\( >0.05 \\) | \\( >0.02 \\) | \\( >0.01 \\) | \\( >0.005 \\) | \\( >0.002 \\) | \\( >0.001 \\) | それ以下 |
|---|---|---|---|---|---|---|---|---|---|---|
| \\( f \\) | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |

- 更新量の成分ごとの最大が 0.2 bohr を超えたら、全体を縮めます（構造最適化の 0.5 bohr と異なる）。
- 端の像は動かしません。`tools/check_geometry_md.py` で、2 次元の模型ポテンシャル上でこの方法が収束し、像の間隔がそろい、最高点の像が鞍点のエネルギーに近づくことを確認しています。

### 8.4 並列化・出力・再開〔公開〕〔規約〕

〔公開〕マニュアル 51.3・51.4・51.7 節：

- `MD.NEB.Parallel.Number` は像を割り当てる MPI プロセスの群の数です。指定しない場合、まず像について並列化し、プロセスが像より多いときに各像の中も並列化します。
  1 にすると像の並列化をせず、全プロセスで像を 1 つずつ計算します（メモリーが最小）。〔規約〕群の数がプロセス数か像の数を超える場合は自動（0）に戻します。
- `System.Name.neb.opt`：最適化の履歴（反復、SD の前因子、最大の力、最大変位、\\( \lVert\mathbf G\rVert \\)、像のエネルギーの和）。
  〔規約〕最大変位の列は見出しが Å ですが、値は bohr です。
- `System.Name.neb.ene`：像ごとの「番号（0 と \\( P+1 \\) が端）、全エネルギー、隣の像との距離、端からの累積距離」。
  〔規約〕2・3 列目は見出しでは距離（bohr）ですが、値は **距離の 2 乗**（bohr²）とその累積です。さらに 5 列目に、原子 2 と原子 1 の x 座標の差（Å）が出ます。
- `System.Name.neb.xyz`：各像の座標（XYZ 形式、見出し行に像の番号とエネルギー）。〔規約〕`System.Name.neb.md` も書かれ、形式は `.md`（10 節）と同じで、力の列は \\( -\mathbf G\_p \\)、速度の列は 0 です。
- `System.Name.dat#`：全像の座標（`NEBp.Atoms.SpeciesAndCoordinates`）と `scf.restart on` を含む再開用の入力。
  〔規約〕両端のエネルギーは `System.Name_0_rst/` に保存され、再開時は両端の SCF を省きます。

## 9. SCF との関係〔規約〕

- ステップ間で引き継ぐもの：前のステップの密度（再開）、密度の外挿の履歴（[密度行列と電荷ミキシング](density-mixing.md) 6・7 節）、
  格子の原点（[実空間グリッド](grid-poisson-xc.md) 1.4 節）。ミキシングの重みは各ステップの最初に `scf.Max.Mixing.Weight` に戻します（同章 4 節）。
- 各ステップで、原子の移動に合わせて原子対の集合（[疎ブロック形式](../design/sparse-block-format.md)）と、
  二中心積分などの行列要素を作り直します。格子の最適化と `DF` では、実空間格子の分割数を最初のステップの値に固定し、
  格子ベクトルの変化に合わせて刻み幅だけを変えます（実効的なカットオフが変わる）。NPT での扱いは〔未確定〕です。
- SCF が収束しなかったステップでも、そのときの力で更新を続けます。`.md` の見出し行の `scf_conv=` が 0 になります。
- 力は SCF の最後の密度行列での解析的な微分です（[全エネルギーと力](energy-force.md)）。MD の保存則の精度は SCF の収束判定に依存するので、MD では `scf.criterion` を十分小さくします。

## 10. 出力ファイル〔公開〕〔規約〕

### `.md`

〔公開〕マニュアル 18.6 節：ステップごとに「原子数」の行、見出し行、原子ごとの行（14 列）のブロックが続きます。
見出し行には自由な文字列の後に `Cell_Vectors=`（\\( \mathbf a\_1,\mathbf a\_2,\mathbf a\_3 \\) の x, y, z の 9 個、Å）と `scf_conv=`（1 で SCF 収束）が付きます。
原子の行は「元素記号、x、y、z（Å）、\\( F\_x, F\_y, F\_z \\)（Hartree/bohr）、\\( v\_x, v\_y, v\_z \\)（m/s）、正味の電荷（電子の電荷を負とする）、磁気モーメント（\\( \mu\_B \\)）、スピンの向きの角度 \\( \theta, \phi \\)（度）」です。

〔規約〕

- 見出し行の自由な部分：MD（`NVE`、`NVT_VS`、`NVT_NH`、`NVT_VS2`、`NVT_Langevin`、`NPT_*`）では
  `time= `（\\( (n-1)\Delta t \\)、fs）、`Energy= `（`Utot`）、`Temperature= `、`(Given Temp.= ` \\( T\_0 \\)`)`。
  `NVT_VS4` では温度の代わりに群ごとの `Temperature[g]= `。それ以外（構造最適化など）は `time=` と `Energy=` だけ。
- `.md` は SCF の直後、座標の更新の **前** に書くので、座標と力は同じ構造に対応します。速度の列は、その時点で保持している速度
  （`NVE` では半ステップの \\( \mathbf v\_{n-1/2} \\)、`NVT_VS` 系では前のステップの更新後の \\( \mathbf u \\)）です。
  構造最適化と `NPT_*` では速度の列は 0 です。温度は 2 節の \\( T\_B \\)（`NPT_*` と `NVT_VS4` を除く）。
- 力の列は、MD の総和 0 の補正の **前** の値です（補正は `.md` の書き出しの後に行う）。
- 正味の電荷と磁気モーメントは原子ごとの電子数（上向き・下向き）から計算されますが、どの段階の値（初期配置か、SCF 後の電荷解析か）かは〔未確定〕です。照合の対象から外します〔設計〕。
- `System.Name.md2` には最後のステップの構造だけが書かれます（上書き）。

### `.ene`

〔規約〕`MD.Type` が `NOMD`、構造最適化、格子の最適化、`EvsLC`、`DF`、`NVE`、`NVT_*` のとき、ステップごとに 1 行を追記します（34 列、空白区切り）。
再実装との互換に必要な列は次のとおりです。それ以外の列は既存実装の内部の量で、照合の対象にしません。

| 列 | 内容 |
|---|---|
| 1 | ステップ \\( n \\) |
| 2 | 時刻 \\( (n-1)\Delta t \\)（fs。構造最適化でも `MD.TimeStep` を使って書く） |
| 14 | `Ukc`（Hartree） |
| 15 | `Utot`（Hartree） |
| 16 | `Utot` + `Ukc` |
| 17 | 化学ポテンシャル（Hartree）。32 列目も同じ |
| 18 | 目標温度 \\( T\_0 \\)（K） |
| 19 | 温度 \\( T\_A \\)（K） |
| 20、21、22 | \\( \zeta \\)、\\( \zeta' \\)、Nosé–Hoover の保存量（`NVT_NH`） |
| 23–25、26 | \\( \mathbf a\_1 \\) の x, y, z と \\( \mathbf a\_3 \\) の z（Å） |
| 30 | 最大の力（そのステップの方式での定義、4.1 節） |
| 31 | 全スピン磁気モーメントの 2 倍 |

〔公開〕マニュアル 20.1 節は 1、2、14–17 列を説明しています。

### その他

- `System.Name.out`：構造最適化の各ステップの「反復、前因子、最大の力、最大変位（Å）、`Utot`」（格子の最適化では `Enpy` と体積（Å³）が加わる。マニュアル 19.4 節）が
  「History of geometry optimization」／「History of cell optimization」の表として、
  最後の構造の座標と力が `<coordinates.forces` ブロックとして含まれます。〔規約〕準 Newton 法の行の前因子の列は最急降下の前因子のままです（意味を持たない）。
- `System.Name.npt.ene`：6.7 節。`System.Name.neb.*`：8.4 節。`System.Name.EvsLC`、`System.Name.DF`：5.5 節。
- `System.Name.dat#`：4.11・6.8・8.4 節。

## 11. 端の場合〔規約〕

- `MD.maxIter` = 1（既定）では、構造最適化・MD を指定しても 1 回の SCF だけで終わります（座標は更新しない）。
- 構造最適化の最初のステップで既に収束していれば、そのまま終わります。
- GDIIS の係数行列が特異（例：同じ勾配が 2 回続いた）だと「収束」として終わります（4.1 節）。
- 曲率条件の破れ：`RF` は分母を \\( 10^{-10} \\) に置き換えて更新を続け、`BFGS` と NEB は無条件に更新し（分母が 0 に近いと発散）、`EF` は更新を飛ばします。
- `Schlegel` の Hessian で \\( r \le B \\) の原子対があると、負の曲率や発散が生じます。
- 熱浴付きの MD を乱数の初速度で始め、固定成分がある場合、固定成分の速度が 0 にならず、`NVT_VS` の温度に入り続けます。
- `NVT_VS2` で全成分が固定された原子は \\( T\_I = 0 \\) でスケーリング係数が定義されません（掛ける成分がないので結果には影響しない）。
- `NVT_NH` で `MD.TempControl` の最初のステップが 1 より大きい場合、\\( n = 1 \\) の予測で使う目標温度は未設定（0 と同じ）です。
- NEB で隣り合う像のエネルギーが等しいと、接線が定義されないことがあります（8.2 節）。
- `MD.CellPressureControl` の範囲を超えた NPT のステップでは外圧が未定義です（6.7 節）。
- 原子がちょうど分数座標 1 にある場合は折り返されません（7 節）。

## 実装への要求〔設計〕

- **構成**：「力・応力を返す計算器」と「状態を進める積分器・最適化器」を分けます。積分器・最適化器は
  \\( (\mathbf x, A, \mathbf g, \Sigma, U) \\) と自分の状態（速度、熱浴、履歴、Hessian）だけを受け取る純粋な関数にし、
  SCF・格子・ファイル入出力に依存させません。これにより、解析的な模型ポテンシャルで単体試験ができます。
- **状態の保存**：最適化器・積分器の状態を 1 つの構造体として直列化し、再開ファイルに書きます（既存実装の `_rst` ディレクトリの形式に合わせる必要はない）。
  再開した計算が、途切れずに計算した場合とステップごとに一致することを試験します。
- **既定の方式と互換モード**：次の点で、既存実装と挙動が異なる「推奨」と、既存実装に合わせる「互換」を切り替えられるようにします。

| 項目 | 推奨（既定） | 互換 |
|---|---|---|
| 速度スケーリングの積分 | 速度 Verlet ＋ スケーリング | 6.4 節の形 |
| 速度の平均の除去 | 全運動量を 0（質量の重み付き） | 単純平均 |
| RF の BFGS で曲率条件が破れたとき | 更新を飛ばす | 分母を \\( 10^{-10} \\) に |
| GDIIS が特異なとき | 履歴を捨てて継続 | 収束扱いで終了 |
| NPT の NH の保存量 | 実際の運動方程式に対応する量 | \\( \tfrac12 Q\zeta^2 \\) を使う |
| WV の \\( \dot F \\) | \\( \dot{\boldsymbol\sigma}\_i\cdot\boldsymbol\sigma\_j + \boldsymbol\sigma\_i\cdot\dot{\boldsymbol\sigma}\_j \\) | \\( 2\dot{\boldsymbol\sigma}\_i\cdot\boldsymbol\sigma\_j \\) |
| 孤立系での折り返し | しない | する |
| `.neb.ene` の距離 | 距離 | 距離の 2 乗 |

- **温度の定義**：\\( T = 2K/(k\_B N\_{\mathrm{dof}}) \\) を一つの関数にし、\\( N\_{\mathrm{dof}} \\) を「\\( 3N \\)」「自由成分の数」「自由成分の数 − 3（重心の拘束）」から選べるようにします。出力ごとにどれを使うかは表で管理します（2 節）。
- **乱数**：Langevin と初速度の乱数は、種を入力で指定できる、環境に依存しない生成器にします。
- **単位**：入力・出力の換算係数は 2 節の値を一か所に定義します。
- **拘束**：固定フラグを、デカルト成分（原子座標の最適化・MD）と分数成分（`RFC*`）の 2 種類として型で区別します。
- **NEB**：像ごとの SCF を独立に呼べるようにし（並列化は外側の責任）、端点の SCF の結果を保存・再利用します。周期系では、端点の座標を最小像で対応付ける前処理を選択肢として持ちます。

## テスト〔設計〕

`tools/check_geometry_md.py` で、下の「解析」「不変量」の項目のうち、模型で書けるものを確認しています。

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 単位換算（fs、u、m/s、GPa）と独立な定数からの値 | 相対 \\( 10^{-5} \\)（GPa は \\( 10^{-4} \\)） |
| 解析 | 調和振動子の速度 Verlet：エネルギー誤差が有界で、\\( \Delta t \\) を半分にすると 1/4 | 比 \\( 4 \pm 0.2 \\) |
| 解析 | 6.4 節の積分の特性根 \\( \lvert z\rvert = \sqrt{1-(\omega\Delta t)^2/2} \\) と、シミュレーションでのエネルギーの減衰率 | 減衰率の相対 5% |
| 解析 | Nosé–Hoover（6.5 節の離散化）：拡張系のエネルギーの変動が \\( O(\Delta t^2) \\) | 比 \\( 4 \pm 1 \\)、変動は \\( 0.05\\,k\_BT \\) 以下（\\( \Delta t = 0.5 \\) fs） |
| 解析 | NPT_NH 型の熱浴：保存量が \\( \tfrac32 Q\zeta^2 \\) を含む（高精度の常微分方程式の解で） | \\( 10^{-6}\\,k\_BT \\) |
| 解析 | 速度スケーリング後の温度が \\( T\_0 + (T-T\_0)\alpha \\) | \\( 10^{-9} \\) K |
| 解析 | GDIIS の係数：縁付き方程式と閉じた形 | \\( 10^{-10} \\) |
| 解析 | BFGS の割線条件 \\( H\_{\mathrm{new}}\Delta\mathbf x = \Delta\mathbf g \\)、\\( H^{-1}\_{\mathrm{new}}\Delta\mathbf g = \Delta\mathbf x \\)、両者が互いに逆行列 | \\( 10^{-9} \\) |
| 解析 | RF の拡大 Hessian の固有ベクトルの step が \\( -(H-\lambda)^{-1}\mathbf g \\)、\\( \lambda = \mathbf g\cdot\Delta\mathbf x < \lambda\_{\min}(H) \\) | \\( 10^{-10} \\) |
| 解析 | 2 次関数（9 次元）で `Opt`、`DIIS`、`BFGS`、`RF`、`EF` がすべて `MD.Opt.criterion` まで収束 | 収束すること |
| 解析 | \\( S = \Omega A^{-\mathsf T} \\)、\\( \partial\Omega/\partial A = S \\) | \\( 10^{-7} \\) |
| 解析 | NEB：2 次元の模型で収束し、像の間隔の相対標準偏差が小さく、最高点が鞍点のエネルギーに近い | 2%、0.05 |
| 不変量 | 初速度の全運動量が 0、温度が \\( T\_1 \\) | \\( 10^{-12} \\) |
| 不変量 | モデル Hessian が対称で、原子対の部分が並進で 0 | \\( 10^{-12} \\) |
| 不変量 | 折り返し：\\( s = 1.25 \to 0.25 \\)、\\( -0.25 \to 0.75 \\)、\\( 1.0 \\) は保持 | 厳密 |
| 不変量 | 構造最適化で、固定成分の座標がビット単位で変わらない | 厳密 |
| 不変量 | 再開：\\( n \\) ステップで止めて `.dat#` から再開した計算と、続けて計算した場合の座標・速度（同じ方式・同じ乱数の種） | \\( 10^{-10} \\) bohr（SCF の初期密度が同じ場合） |
| 不変量 | NVE の実 DFT 計算（H₂O 分子、\\( \Delta t = 0.5 \\) fs、`scf.criterion` \\( 10^{-10} \\)）で `Utot` + `Ukc` の変動 | 100 ステップで \\( 10^{-4} \\) Hartree 程度（傾向の確認） |
| 照合 | 既存実装の fixture：同じ初期構造・同じ力の列（力を記録したもの）を与えたとき、`Opt`、`DIIS`、`BFGS`、`RF`、`EF` の 1 ステップの更新量 | \\( 10^{-8} \\) bohr |
| 照合 | 同様に、`NVE`、`NVT_VS`（互換モード）、`NVT_NH` の 1 ステップの座標・速度・\\( \zeta \\) | \\( 10^{-10} \\)（相対） |
| 照合 | `OptC1`、`OptC5`、`RFC5` の 1 ステップの格子と座標（応力と力を与えて） | \\( 10^{-8} \\) bohr |
| 照合 | `Schlegel` の初期 Hessian（H₂O、Si 結晶） | \\( 10^{-10} \\) |
| 回帰 | H₂O 分子の `RF` 最適化：収束した構造（結合長・角度）と反復回数 | 結合長 \\( 10^{-3} \\) Å、角度 0.1°、反復回数は目安 |
| 回帰 | Si 結晶の `RFC5`：収束した格子定数 | \\( 10^{-3} \\) Å |
| 回帰 | NEB：マニュアルの例（C₂H₄ の回転など）で障壁の高さ | 10 meV |

照合のための fixture は、既存実装に「力を外から与える」経路がないため、SCF を含む小さな系で 2–3 ステップ実行し、
各ステップの座標・力・応力（`.md`、`.out`）を記録したものから、1 ステップの更新を再現する形で作ります。

## 〔未確定〕

- `MD.Opt.Init.Hessian FF`（力場の当てはめによる初期 Hessian）の式。公開資料に記載がなく、使用例も見当たらないため、再実装の対象外とする案。
- `RFC5`–`RFC7` の `Schlegel` 初期 Hessian の、格子どうし・原子と格子の間の部分の正確な式（収束の経路だけに影響）。
- NEB での `MD.Opt.DIIS.History` の実効的な既定値。NEB の入力の読み込みは 3 を既定にしますが、各像の SCF の入力処理で 2 に上書きされるように見え、実行して確認する必要があります。
- WV の `Ubaro`（`.npt.ene` の 8 列目）の定義と、`NPT.WV.F0` の用途（運動方程式には使われていない）。
- `.ene` の 3–13 列、27–29 列、33–34 列の意味（内部の量。互換の対象にするか）。
- NPT で `MD.Fixed.XYZ` を与えた場合の挙動（マニュアルは併用不可とするが、入力の検査はない）。
- 本家 OpenMX（CPU 版）との差：参照したのは GPU 版の派生で、`NVT_VS` 系の積分の形（6.4 節）や既定値が本家でも同じかは未確認です。
- `Langevin.Friction.Factor` の単位（原子単位の時間の逆数と判断したが、公開資料の裏付けがない）。
