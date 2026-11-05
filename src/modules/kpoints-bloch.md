# k 点と Bloch 和

> 状態：試作（レビュー待ち）。Bloch 和の式は技術ノート 1（Eq. 1、4、28–29）をこのノートの記法で書き直したもの、
> k 点の生成と重みは既存実装の挙動として確認したものです。
> 数値的な性質は `tools/check_kpoints_eigen.py` で確認しています。

## この章の範囲

実空間の行列 \\( O(\mathbf L) \\)（[疎ブロック形式](../design/sparse-block-format.md)）から、
各 k 点の行列 \\( H(\mathbf k) \\)、\\( S(\mathbf k) \\) を作る規約と、k 点の集合と重みを定めます。
k 点ごとの固有値問題の解き方は [一般化固有値問題と S の正則化](generalized-eigen.md)、
固有ベクトルから密度行列を作る詳細は [密度行列と電荷ミキシング](density-mixing.md) で扱います。

## 格子と逆格子〔設計〕

単位胞のベクトルを \\( \mathbf a\_1, \mathbf a\_2, \mathbf a\_3 \\)、逆格子ベクトルを
\\( \mathbf b\_i\cdot\mathbf a\_j = 2\pi\delta\_{ij} \\) で定めます。
格子ベクトルは整数の 3 つ組 \\( \mathbf n = (n\_1, n\_2, n\_3) \\) で
\\( \mathbf L = \sum\_i n\_i \mathbf a\_i \\)、k 点は分数座標 \\( \mathbf q = (q\_1, q\_2, q\_3) \\) で
\\( \mathbf k = \sum\_i q\_i \mathbf b\_i \\) と表します。このとき

$$
\mathbf k\cdot\mathbf L = 2\pi\\,(q\_1 n\_1 + q\_2 n\_2 + q\_3 n\_3)
$$

です。k 点の入出力と内部表現は分数座標 \\( \mathbf q \\) を基本とします。

## Bloch 和〔公開〕〔規約〕

〔公開〕技術ノート 1 の Eq. 1 は、波数 \\( \mathbf k \\)、バンド \\( \mu \\) の Kohn–Sham 軌道を

$$
\psi\_{\mu\mathbf k}(\mathbf r) = \frac{1}{\sqrt{N\_{\mathrm{cell}}}}\sum\_{\mathbf L} e^{i\mathbf k\cdot\mathbf L}
\sum\_{I\alpha} c\_{I\alpha,\mu}(\mathbf k)\\, \phi\_{I\alpha}(\mathbf r - \mathbf R\_I - \mathbf L)
$$

と展開します（規格化の係数の書き方は原文と異なりますが、\\( c \\) の規格化 \\( c^\dagger S c = 1 \\) と整合するように選んであります）。
位相 \\( e^{i\mathbf k\cdot\mathbf L} \\) には **格子ベクトルだけ** が入り、原子位置 \\( \mathbf R\_I \\) は入りません。
これを Hamiltonian に代入すると（同 Eq. 28–29）、

$$
O\_{I\alpha, J\beta}(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L}\\, O\_{I\alpha, J\beta}(\mathbf L),\qquad
O\_{I\alpha, J\beta}(\mathbf L) = \langle \phi\_{I\alpha}(\mathbf r - \mathbf R\_I) | \hat O | \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L) \rangle
$$

となります。\\( \hat O \\) は重なり、運動エネルギー、Kohn–Sham Hamiltonian などです。
\\( \mathbf L \\) は **列側**（ket 側）の原子 \\( J \\) の周期像を指し、符号は \\( +i \\) です。

〔規約〕既存実装の \\( O(\mathbf k) \\) の \\( (I, J) \\) ブロックは、\\( J \\) の周期像（位置 \\( \mathbf R\_J + \mathbf L \\)）に
わたる \\( e^{2\pi i\\,\mathbf q\cdot\mathbf n}\\, O\_{IJ}(\mathbf L) \\) の和です（複数の周期像は加算されます）。
つまり上の式どおりで、原子位置は位相に入りません。

### 性質〔公開〕

以下は定義から直ちに従い、テストの不変量として使えます。

- Hermite 性：実空間で \\( O\_{J\beta, I\alpha}(-\mathbf L) = \overline{O\_{I\alpha, J\beta}(\mathbf L)} \\) が成り立てば、\\( O(\mathbf k) \\) は Hermite です。
- 周期性：位相が格子ベクトルだけなので、逆格子ベクトル \\( \mathbf G \\) について \\( O(\mathbf k + \mathbf G) = O(\mathbf k) \\) です。
  k 点の分数座標を 1 ずらしても行列は変わりません。
- 時間反転：\\( O(\mathbf L) \\) が実行列なら（共線計算の \\( H \\) と \\( S \\)）、\\( O(-\mathbf k) = \overline{O(\mathbf k)} \\) です。
  したがって固有値は \\( \varepsilon\_\mu(-\mathbf k) = \varepsilon\_\mu(\mathbf k) \\)、固有ベクトルは \\( c(-\mathbf k) = \overline{c(\mathbf k)} \\) に選べます。

### 原子位置を位相に含める規約との関係〔公開〕

位相を \\( e^{i\mathbf k\cdot(\mathbf R\_J + \mathbf L - \mathbf R\_I)} \\) とする規約（原子位置を含める規約）では、
対角行列 \\( P(\mathbf k) = \mathrm{diag}\left(e^{i\mathbf k\cdot\mathbf R\_I}\right) \\)（原子 \\( I \\) の軌道すべてに同じ位相）を使って

$$
\tilde O(\mathbf k) = P(\mathbf k)^\dagger\\, O(\mathbf k)\\, P(\mathbf k),\qquad \tilde c(\mathbf k) = P(\mathbf k)^\dagger c(\mathbf k)
$$

となります。固有値は変わりません。固有ベクトルの成分は原子ごとの位相だけ変わります。
\\( \tilde O(\mathbf k + \mathbf G) \\) は \\( \tilde O(\mathbf k) \\) と一致しない（ユニタリー同値）ので、
k による微分（速度演算子、Berry 位相など）を扱う後処理では、どちらの規約かを必ず明示する必要があります。
このノートでは既存実装と同じく、格子ベクトルだけの規約を基本とします。

### 実空間の密度行列〔公開〕〔規約〕

〔公開〕技術ノート 1 の Eq. 4 に対応して、k 点の固有ベクトルから実空間の密度行列を

$$
\rho\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_\mu f\_{\mu\mathbf k}\\,
e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\mu}(\mathbf k)}\\, c\_{J\beta,\mu}(\mathbf k)
$$

と作ります（\\( \sum\_{\mathbf k} w\_{\mathbf k} = 1 \\)、\\( f \\) は占有数）。この位相の向きで、
バンドエネルギーは \\( \sum\_{\mathbf L}\sum\_{I\alpha J\beta} \rho\_{I\alpha,J\beta}(\mathbf L)\\, H\_{I\alpha,J\beta}(\mathbf L) \\)、
電子数は同じ式の \\( H \\) を \\( S \\) に替えたものになります。

〔規約〕共線計算の既存実装は、上の式の実部だけを足し込みます。時間反転で \\( \mathbf k \\) と \\( -\mathbf k \\) の寄与は互いに複素共役なので、
\\( \pm\mathbf k \\) の対を重み 2 の 1 点で代表させても実部をとれば厳密です（下の「k 点の削減」）。

## 計算の種類と行列の型〔公開〕〔規約〕

〔公開〕マニュアル（キーワードの節）によれば、`scf.EigenvalueSolver` で解法を選び、
`Cluster`（孤立系）、`Band`（周期系、k 点サンプリング）、`DC`、`DC-LNO`、`Krylov`（O(N) 法）、低次スケーリング法などがあります。
〔規約〕マニュアル 8.2 節は低次スケーリング法の値を `ON2` と書いていますが、既存実装が受け付ける値は `Cluster2` です（[O(N) 法](on-methods.md)）。
`scf.SpinPolarization` は `Off`（非スピン偏極）、`On`（共線スピン偏極）、`NC`（非共線）を取り、
`scf.SpinOrbit.Coupling` は非共線計算でのみ使えます。

〔規約〕既存実装で各 k 点（または各局所問題）が解く行列の型は次のとおりです。
\\( N \\) は単位胞（または計算領域）の全軌道数です。

| 解法 | スピン | 解く行列 | 型 | 個数 |
|---|---|---|---|---|
| Cluster | `Off` / `On` | \\( H(\Gamma) \\)、\\( S(\Gamma) \\) | 実対称 \\( N\times N \\) | スピンごとに 1 つ（`Off` は 1 つ） |
| Cluster | `NC` | \\( H(\Gamma) \\)、\\( S(\Gamma) \\) のスピノル版 | 複素 Hermite \\( 2N\times 2N \\) | 1 つ |
| Band | `Off` / `On` | \\( H(\mathbf k) \\)、\\( S(\mathbf k) \\) | 複素 Hermite \\( N\times N \\) | k 点 × スピン |
| Band | `NC` | スピノル版 | 複素 Hermite \\( 2N\times 2N \\) | k 点ごとに 1 つ |
| DC | `Off` / `On` | 原子ごとの局所行列 | 実対称 | 原子 × スピン |
| DC | `NC` | 原子ごとの局所行列のスピノル版 | 複素 Hermite | 原子ごと |

- Cluster の \\( H(\Gamma) \\) は \\( \sum\_{\mathbf L} O(\mathbf L) \\)（\\( \mathbf k = 0 \\) の Bloch 和）です。
  孤立分子のように周期像の近接がなければ \\( \mathbf L = 0 \\) だけが残ります。周期的な単位胞で Cluster を指定した場合は、
  周期像の寄与がそのまま加算されて Γ 点の計算になります。
- `Band` で `scf.Kgrid` が `1 1 1` のとき、既存実装は解法を Cluster（周期境界条件つきの Γ 点計算）に切り替え、実数演算で解きます。
- 共線スピン偏極では、上向き・下向きスピンの \\( H \\) は別の行列で、\\( S \\) は共通です。
- 共線計算の \\( O(\mathbf L) \\) はすべて実数です。非共線計算では、スピン非対角ブロックと、
  スピン軌道相互作用・LDA+U・スピン拘束・Zeeman 項などによる虚部が加わります。

〔規約〕O(N) 解法（`DC`、`Krylov`、`DC-LNO`、`Cluster2`）の行列の型と扱いは [O(N) 法](on-methods.md)、`NEGF` の行列（中心領域と電極の k 空間ブロック）は [NEGF による電子輸送](negf.md) で定めます。
〔未確定〕`Cluster-LNO` の行列の型。
既存実装には `Recursion`、`Cluster-DIIS`、`Cluster2`、`EGAC` という値も残っていますが、マニュアルの記載とあわせて別途確認が必要です。

## 非共線スピノルの並び〔規約〕

非共線計算の \\( 2N \\) 次元の行列は **スピンが外側、軌道が内側** の並びです。行・列の番号を

$$
p = \sigma N + a,\qquad \sigma \in \\{0\\,(\uparrow),\\, 1\\,(\downarrow)\\},\quad a = 0, \dots, N-1
$$

とします。\\( a \\) は共線計算と同じ全軌道の番号（原子を入力順に並べ、各原子の中は [PAO 基底](pao-basis.md) の順）です。
したがって

$$
H^{(2)}(\mathbf k) = \begin{pmatrix} H^{\uparrow\uparrow}(\mathbf k) & H^{\uparrow\downarrow}(\mathbf k) \\\\ H^{\downarrow\uparrow}(\mathbf k) & H^{\downarrow\downarrow}(\mathbf k) \end{pmatrix},\qquad
S^{(2)}(\mathbf k) = \begin{pmatrix} S(\mathbf k) & 0 \\\\ 0 & S(\mathbf k) \end{pmatrix} = I\_2 \otimes S(\mathbf k)
$$

です。スピンの量子化軸は全体座標の \\( z \\) 軸で、\\( H^{\uparrow\downarrow}\_{ab} = \langle a\uparrow | \hat H | b\downarrow\rangle \\) です。

- 各スピンブロックは、それぞれの実空間行列 \\( H^{\sigma\sigma'}(\mathbf L) \\) から同じ Bloch 和で作ります。
- 既存実装は \\( \uparrow\uparrow \\)、\\( \downarrow\downarrow \\)、\\( \uparrow\downarrow \\) の 3 ブロックを Bloch 和で作り、
  \\( H^{\downarrow\uparrow}(\mathbf k) = H^{\uparrow\downarrow}(\mathbf k)^\dagger \\) で残りを埋めます。
- 実空間では \\( H^{\sigma'\sigma}\_{J\beta, I\alpha}(-\mathbf L) = \overline{H^{\sigma\sigma'}\_{I\alpha, J\beta}(\mathbf L)} \\) が成り立ちます。
- 重なりはスピンに依存しないので、\\( S^{(2)} \\) は \\( S \\) を 2 回並べたブロック対角です。
- 固有ベクトル \\( c \\)（長さ \\( 2N \\)）は、前半 \\( N \\) 成分が上向きスピン、後半 \\( N \\) 成分が下向きスピンの係数です。

〔規約〕非共線計算の実空間密度行列（スピン 2×2 成分）は、上の Bloch 和の規約から自然に導かれる
\\( \rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L) = \sum w f e^{i\mathbf k\cdot\mathbf L}\\, \overline{c\_{\sigma I\alpha}}\\, c\_{\sigma' J\beta} \\) の形です。
成分の保持方法は [密度行列と電荷ミキシング](density-mixing.md) を参照してください。

## k 点メッシュ〔公開〕〔規約〕

〔公開〕マニュアル 4.0（15 章）によれば、`scf.Kgrid` で逆格子ベクトル方向の分割数 \\( (N\_1, N\_2, N\_3) \\) を、
`scf.Generation.Kpoint` で生成方法を指定します。選べるのは `GCenter`（既定）と `Regular` です。
`GCenter` は分割数の偶奇によらず Γ 点を含み、\\( \mathbf k \\) と \\( -\mathbf k \\) を対等に扱い、
\\( \mathbf k \equiv -\mathbf k \pmod{\mathbf G} \\) を満たす時間反転不変運動量（TRIM）を含みます。

〔規約〕`scf.Kgrid` を省略したときの既定値は `4 4 4` です（マニュアルは Band で必須としています）。
キーワードの値は大文字・小文字を区別しません。

### 座標

各方向 \\( i \\)、番号 \\( l = 0, \dots, N\_i - 1 \\) について、分数座標は次のとおりです。

| 方式 | \\( q\_{i,l} \\) | 備考 |
|---|---|---|
| `GCenter` | \\( l / N\_i \\) | \\( [0, 1) \\) に並ぶ。\\( l = 0 \\) が Γ。\\( N\_i \\) が偶数なら \\( q = 1/2 \\)（境界）も含む |
| `Regular` | \\( -\tfrac12 + \tfrac{2l + 1}{2N\_i} \\) | Γ を含むのはすべての \\( N\_i \\) が奇数のときだけ。境界上の点を含まない |

- どちらの方式でも \\( N\_i = 1 \\) の方向は \\( q\_i = 0 \\) です（`Regular` の式もこの場合 0 になります）。
- 〔公開〕`Regular` の式は Monkhorst–Pack の式（PRB **13**, 5188 (1976)）の
  \\( u\_r = (2r - N - 1)/(2N) \\)、\\( r = 1, \dots, N \\) と、\\( r = l + 1 \\) として同一です。
  したがって `Regular` は「Monkhorst–Pack 型の格子（記号の縮約はしない）」に当たります。
- `GCenter` の点は \\( [0, 1) \\) にありますが、周期性 \\( O(\mathbf k + \mathbf G) = O(\mathbf k) \\) により、
  \\( (-\tfrac12, \tfrac12] \\) に折り返した Γ 中心の格子と同じ行列を与えます。
  これは \\( N\_1\times N\_2\times N\_3 \\) 超格子の Γ 点計算と同じ固有値の集合を与えます（テスト参照）。
- 〔規約〕既存実装は、分割数が 2 以上の方向の分数座標に、方向 1・2・3 の順に
  \\( +\delta \\)、\\( -\delta \\)、\\( +2\delta \\)（\\( \delta = 10^{-12} \\)）を加えます。固有値ルーチンの安定化のための微小なずれで、
  再実装で再現する必要はありません（[ビット一致を求めない](../process/writing-rules.md)）。

### k 点の削減と重み〔規約〕

共線計算（`scf.SpinPolarization` が `Off` または `On`）の `Band` では、時間反転による
\\( \varepsilon(\mathbf k) = \varepsilon(-\mathbf k) \\) だけを使って k 点を減らします。結晶の点群対称性は使いません。

- 格子点 \\( \mathbf l = (l\_1, l\_2, l\_3) \\) の相手 \\( \bar{\mathbf l} \\) は、\\( -\mathbf q \\) に（\\( \mathbf G \\) を除いて）対応する格子点です。
  - `Regular`：\\( \bar l\_i = N\_i - 1 - l\_i \\)
  - `GCenter`：\\( \bar l\_i = (N\_i - l\_i) \bmod N\_i \\)（\\( l\_i = 0 \\) と \\( 2l\_i = N\_i \\) は自分自身が相手）
- \\( \bar{\mathbf l} = \mathbf l \\) の点（TRIM）は重み 1、それ以外は対の一方を重み 2 で残し、他方を除きます。
  残すのは、\\( l\_1 \\) を最も外側、\\( l\_3 \\) を最も内側とする辞書式の順で先に現れる方です。
- 重みの総和は \\( N\_1 N\_2 N\_3 \\) です。密度行列などで使う重みは、これで割って総和 1 に規格化します。
- TRIM の数：`GCenter` では各方向で \\( N\_i \\) が偶数なら 2 個、奇数なら 1 個の積です。
  `Regular` では、すべての \\( N\_i \\) が奇数のときの Γ だけ（1 個）、そうでなければ 0 個です。

非共線計算の `Band` では削減をせず、\\( N\_1 N\_2 N\_3 \\) 個の点をすべて重み 1 で使います。
非共線の \\( H(\mathbf L) \\) は複素数を含むので \\( H(-\mathbf k) = \overline{H(\mathbf k)} \\) が成り立たず、
時間反転にはスピン反転を伴う反ユニタリー変換が必要で、磁性があればそれも破れるためです〔公開〕（一般的な事実）。

〔未確定〕

- 既存実装には、原子配置の対称性を使って Monkhorst–Pack 格子を縮約する別の方式（許容誤差のキーワード `scf.MP.criterion`）が残っています。
  マニュアル 4.0 は `GCenter` と `Regular` しか記載しておらず、既存実装でこの方式を選ぶ値もマニュアルにありません。
  共線の `Band` でのみ使え、NEGF と非共線では使えません。このノートでは対象外とします。
- バンド図（k 点の経路）の生成規約（[解析機能](analysis.md) 2 節）。

〔規約〕SCF 以外で使う k 点は [解析機能](analysis.md) で定められています。状態密度（`Dos.Kgrid`）は \\( q\_{i,l} = l/N\_i \\) の Γ 中心の格子で、`scf.Generation.Kpoint` によらず、時間反転などによる削減はしません。
交換結合定数（jx）は Monkhorst–Pack 型の格子で削減しません。光学伝導度（`CDDF.Kgrid`）は Monkhorst–Pack 型で、\\( \pm\mathbf k \\) を重み 2 の 1 点にまとめます。

## 実装への要求〔設計〕

- Bloch 和は、疎ブロック形式の \\( (I, J, \mathbf n) \\) ごとに位相 \\( e^{2\pi i\\,\mathbf q\cdot\mathbf n} \\) を掛けて足し込む、
  状態を持たない関数にします。入力は分数座標 \\( \mathbf q \\)、出力は密行列（列優先）です。
- 位相の規約（格子ベクトルのみ／原子位置を含む）は型またはオプションで明示し、既定は格子ベクトルのみとします。
- 共線計算では \\( O(\mathbf L) \\) が実数であることを使い、Γ 点では実数の行列を返す経路を用意します。
- 非共線計算のスピノル行列は \\( p = \sigma N + a \\) の並びとし、\\( S^{(2)} \\) は明示的に作らず \\( I\_2\otimes S \\) として扱えるようにします。
- k 点の集合は「分数座標と規格化された重みの列」として生成し、生成方法（`GCenter`／`Regular`）と
  時間反転による削減の有無を引数にします。削減は共線計算でのみ許します。
- 既存実装との照合用に、微小なずれ \\( \delta \\) を加えるかどうかをオプションにしても構いません（既定は加えない）。

## テスト

`tools/check_kpoints_eigen.py` は、ランダムな実数の \\( H(\mathbf L) \\)、\\( S(\mathbf L) \\)
（\\( O(-\mathbf L) = O(\mathbf L)^{\mathsf T} \\)）を持つ小さな模型で、以下の多くを確認しています。

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 不変量 | \\( H(\mathbf k) \\)、\\( S(\mathbf k) \\) が Hermite | 丸め誤差 |
| 不変量 | \\( O(\mathbf k + \mathbf G) = O(\mathbf k) \\)（分数座標を整数だけずらす） | \\( 10^{-12} \\) |
| 不変量 | 共線：\\( O(-\mathbf k) = \overline{O(\mathbf k)} \\) | 丸め誤差 |
| 不変量 | 原子位置を含む規約に変えても固有値が不変 | \\( 10^{-11} \\) |
| 不変量 | 削減した k 点と重みでの和が、全格子点での和と一致（\\( f(\mathbf k) = f(-\mathbf k) \\) を満たす量、例：占有状態の固有値和） | \\( 10^{-10} \\) |
| 不変量 | 重みの総和が \\( N\_1 N\_2 N\_3 \\)、TRIM の数が上の規則どおり | 厳密 |
| 不変量 | `GCenter` の \\( N\_1\times N\_2\times N\_3 \\) 格子の固有値の集合が、同じ大きさの超格子の Γ 点の固有値と一致 | \\( 10^{-10} \\) |
| 不変量 | 削減した k 点から作った \\( \rho(\mathbf L) \\) で、\\( \sum \rho H \\) がバンドエネルギーに、\\( \sum \rho S \\) が電子数に一致 | \\( 10^{-11} \\) |
| 不変量 | 非共線：\\( H^{(2)} \\) が Hermite、\\( S^{(2)} = I\_2 \otimes S \\) | 丸め誤差 |
| 不変量 | 非共線でスピン非対角項が 0 のとき、\\( 2N \\) 次元の固有値が共線の 2 スピン分の固有値の和集合になる | \\( 10^{-12} \\) |
| 照合 | k 点の分数座標と重みの一覧を、既存実装の出力（`scf.Kgrid` が `2 3 4`、`3 3 3`、`4 4 1` など、両方式）と比較 | \\( 10^{-10} \\)（\\( \delta \\) を除く） |
| 照合 | 既存実装から書き出した \\( H(\mathbf L) \\)、\\( S(\mathbf L) \\) と、いくつかの k 点の \\( H(\mathbf k) \\)、\\( S(\mathbf k) \\)（共線・非共線） | 要決定 |
