# 密度行列と電荷ミキシング

> 状態：試作（レビュー待ち）。占有数・密度行列・ミキシング・収束判定・外挿の式は、技術ノート 1・8・12 と
> マニュアル 4.0（8.2 節、16 章、17 章）を、このノートの記法で書き直したものです。
> 公開資料に書かれていない規則（既定値、残差ノルムの正規化、重みの更新則、収束判定の組み合わせなど）は
> 既存実装の挙動として確認したもので、〔規約〕と記します。
> 数値的な性質は `tools/check_density_mixing.py` で確認しています。

## この章の範囲

SCF ループのうち、固有値問題を解いた後から次の反復の Hamiltonian を作る前までを扱います。

1. 固有値から占有数と化学ポテンシャルを決める。
2. 固有ベクトルと占有数から実空間の密度行列 \\( \rho(\mathbf L) \\) とエネルギー重み付き密度行列 \\( E(\mathbf L) \\) を作る。
3. \\( \rho \\) から格子上の電子密度 \\( n(\mathbf r) \\) を作る。
4. 入力と出力の差（残差）を使って、次の反復の入力（密度行列、密度、ポテンシャル、または Hamiltonian）を作る（ミキシング）。
5. 収束を判定する。

固有値問題そのものは [一般化固有値問題と S の正則化](generalized-eigen.md)、k 点と Bloch 和は [k 点と Bloch 和](kpoints-bloch.md)、
格子上の積分と Poisson 方程式は [実空間グリッド・Poisson・交換相関](grid-poisson-xc.md) にあります。

## SCF ループの全体像〔公開〕〔規約〕

〔公開〕技術ノート 12 とマニュアル 16.1 節によると、SCF は「入力の密度（または密度行列）から Hamiltonian を作る → 固有値問題を解く →
出力の密度行列を作る → 入力と出力を混ぜて次の入力を作る」の繰り返しで、混ぜる量と方法は `scf.Mixing.Type` で選びます。

〔規約〕既存実装の 1 反復の順序は次のとおりです。

1. 入力の密度 \\( n^{\mathrm{in}} \\) から Hartree・交換相関ポテンシャルを作り、Hamiltonian \\( H(\mathbf L) \\) を組み立てる。
   方式 RMM-DIISV ではこの段階でポテンシャルを、方式 RMM-DIISH ではこの直後に Hamiltonian を混ぜる。
2. 各 k 点・スピンで固有値問題を解く。
3. 化学ポテンシャルを決め、占有数、バンドエネルギー \\( U\_{\mathrm{ele}} \\)、出力の密度行列 \\( \rho^{\mathrm{out}} \\)、エネルギー重み付き密度行列を作る。
4. 収束を判定する（後述）。収束していればループを抜ける。
5. 方式に応じて密度行列または密度を混ぜ、次の反復の入力の密度 \\( n^{\mathrm{in}} \\) を格子上に作る。

収束判定が混合より前にあるので、判定に使う残差ノルムは、直近に混合を行ったときの値（1 つ前の反復の残差）です。

## 1. 占有数と化学ポテンシャル

### 電子数〔公開〕〔規約〕

〔公開〕マニュアル（8.2 節 `scf.system.charge`、31 章）によると、`scf.system.charge` は系の余分な電荷で、正の値は正孔ドープ、負の値は電子ドープ、既定値は 0 です。

〔規約〕占有させる電子数は

$$
N\_{\mathrm{e}} = \sum\_{I} Z\_I^{\mathrm{val}} - Q,\qquad Q = \texttt{scf.system.charge}
$$

です。\\( Z\_I^{\mathrm{val}} \\) は原子 \\( I \\) の擬ポテンシャルの価電子数（芯電荷）です。

### Fermi–Dirac 分布〔公開〕〔規約〕

〔公開〕マニュアル 8.2 節によると、`scf.ElectronicTemperature` は電子温度（K）で、既定値は 300 K です。

〔規約〕占有数は Fermi–Dirac 分布

$$
f(\varepsilon) = \frac{1}{1 + e^{x}},\qquad x = \frac{\varepsilon - \mu}{k\_{\mathrm B} T}
$$

です。他のスメアリング関数（Gauss、Methfessel–Paxton など）はありません。

- \\( k\_{\mathrm B} = 8.617251324\times 10^{-5} \\) eV/K、1 Hartree = 27.2113845 eV で、温度を Hartree 単位に直します。
- 指数の溢れを避けるため、\\( x \\) を \\( [-60, 60] \\) に切り詰めてから評価します（Band の場合）。
- 占有数が \\( 10^{-13} \\) 程度より小さい状態は、密度行列の和から省いてかまいません（結果への影響は丸め誤差程度）。

### 化学ポテンシャル〔規約〕

\\( \mu \\) は電子数の条件

$$
N(\mu) = g\_s \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\sigma}\sum\_{\nu=1}^{N\_{\mathrm{st}}} f\left(\varepsilon\_{\nu\sigma}(\mathbf k)\right) = N\_{\mathrm e}
$$

から決めます。\\( w\_{\mathbf k} \\) は \\( \sum\_{\mathbf k} w\_{\mathbf k} = 1 \\) に規格化した k 点の重み（[k 点と Bloch 和](kpoints-bloch.md)）、
\\( g\_s \\) は非スピン偏極で 2、共線スピン偏極と非共線で 1 です。

- **解法は二分法**です。区間の初期値は \\( [-20, 20] \\) Hartree（Cluster では \\( [-30, 30] \\)）、
  停止条件は \\( \lvert N(\mu) - N\_{\mathrm e}\rvert < 10^{-12} \\)、反復の上限は 1000 回です。
  \\( N(\mu) \\) は単調増加なので、二分法は必ず収束します。
- 〔規約〕この区間と許容誤差は Cluster と Band の値です。O(N) 解法（`DC`、`Krylov`、`DC-LNO`）は解法ごとに別の初期区間・許容誤差を使います（[O(N) 法](on-methods.md) の「電子数と化学ポテンシャル」）。
  NEGF の中心領域では電子数の拘束による探索をせず、電極の化学ポテンシャルから電子数が決まります（Band の段階だけ別の規則。[NEGF による電子輸送](negf.md)）。
  また O(N) 解法の \\( \rho \\) は下の Hermite 性を満たしません（\\( \rho\_{I\alpha,J\beta}(\mathbf L) \ne \rho\_{J\beta,I\alpha}(-\mathbf L) \\)。[O(N) 法](on-methods.md)）。
- 既存実装は、まず 1200 K で \\( \mu \\) を求め、それを最初の試行点にして、指定温度でもう一度区間 \\( [-20, 20] \\) から二分法を行います。
  最終結果は指定温度での電子数の条件だけで決まるので、再実装では 1 回の二分法（または Newton 法との併用）で十分です。
- 共線スピン偏極では、上向き・下向きスピンに **共通の** \\( \mu \\) を使います。スピンごとの電子数は固定せず、全電子数だけを合わせます。
- \\( N\_{\mathrm{st}} \\) は計算する固有状態の数で、既存実装は \\( N\_{\mathrm e}/2 \\) に余分の状態を加えた数（基底の数が上限）を求めます。
  高温で \\( N\_{\mathrm{st}} \\) 番目の状態の占有数が無視できない場合は、状態数を増やす必要があります（〔設計〕：占有数が閾値を超えたら警告）。

> **注意**：絶縁体で \\( k\_{\mathrm B}T \\) がギャップより十分小さいと、ギャップの中で \\( N(\mu) \\) はほぼ平坦です。
> 停止条件は電子数で判定するので、得られる \\( \mu \\) はギャップ内のどこかに止まり、中央とは限りません
> （`tools/check_density_mixing.py`：300 K、ギャップ 0.4 Hartree の 2 準位系で、中央から 0.1 Hartree ずれた点で停止）。
> 占有数と密度は影響を受けませんが、\\( \mu \\) を出力や後処理に使うときは、この不定性を許容誤差に含めます。

### 固定スピンモーメント〔規約〕

既存実装には、共線計算で上向き・下向きスピンの電子数を別々に固定する（スピンごとに別の化学ポテンシャルを使う）入力キーワードは見当たりません。
スピン配置の制御には、非共線計算のスピン方向の拘束（`scf.Constraint.NC.Spin`）や Zeeman 項（`scf.NC.Zeeman.Spin`）を使います。
これらは Hamiltonian に項を加える方法で、占有数の決め方は上と同じです。

### バンドエネルギー〔規約〕

$$
U\_{\mathrm{ele}} = g\_s \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\sigma\nu} f\left(\varepsilon\_{\nu\sigma}(\mathbf k)\right) \varepsilon\_{\nu\sigma}(\mathbf k)
$$

を毎反復計算します。収束判定と、ミキシングの重みの調整に使います。

### 電子温度の自動制御〔規約〕

`scf.Mixing.Control.Temp`（論理値、既定は off）を on にすると、SCF が振動したときに電子温度を一時的に上げます。

- 残差ノルムを \\( \lVert R\rVert \\)（後述の NormRD）とし、直近 5 回の \\( \lVert R\rVert^2 \\) の平均が現在の \\( \lVert R\rVert^2 \\) より小さく（つまり残差が増加傾向で）、
  \\( \lVert R\rVert > 0.01 \\)、反復が 6 回を超えていれば、温度を指定値の 10 倍にします（1 回の SCF で 1 度だけ）。
- その後 \\( \lVert R\rVert \le 0.01 \\) になったら、温度を指定値に戻します。
- 温度を変えるたびに、Pulay 型ミキシングの履歴をリセットして、`scf.Mixing.StartPulay` からの手順をやり直します。

## 2. 実空間の密度行列

### 定義と位相〔公開〕〔規約〕

〔公開〕技術ノート 1 の Eq. 4 と [k 点と Bloch 和](kpoints-bloch.md) の規約
\\( O(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L} O(\mathbf L) \\) に合わせて、密度行列を

$$
\rho^{\sigma}\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\nu} f\_{\nu\sigma}(\mathbf k)\\,
\mathrm{Re}\left[ e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\nu\sigma}(\mathbf k)}\\, c\_{J\beta,\nu\sigma}(\mathbf k) \right]
$$

と作ります。\\( c \\) は \\( c^\dagger S(\mathbf k) c = 1 \\) に規格化した固有ベクトル、\\( \mathbf L \\) は列側（\\( J \\)）の原子の周期像です。

〔規約〕

- 既存実装はこの式のとおり、\\( \cos(\mathbf k\cdot\mathbf L) \\) と \\( \sin(\mathbf k\cdot\mathbf L) \\) を掛けた実部だけを足し込み、
  最後に k 点の重みの総和で割ります（整数の重みを使うので、規格化は最後にまとめて行います）。
- 原子位置は位相に入りません。
- \\( \rho \\) は [疎ブロック形式](../design/sparse-block-format.md) と同じ原子対 \\( (I, J, \mathbf L) \\) の集合の上でだけ作ります。
  \\( (I, J, \mathbf L) \\) と \\( (J, I, -\mathbf L) \\) の両方を持ちます。
- 非スピン偏極では、1 スピン分の密度行列（占有数 \\( f \le 1 \\)）だけを持ちます。全電子数・全密度はその 2 倍です。
- Cluster（Γ 点のみ）では \\( w = 1 \\)、位相は 1 です。

〔設計〕時間反転対称（共線計算で \\( O(\mathbf L) \\) が実数）のとき、\\( \pm\mathbf k \\) の寄与は互いに複素共役なので、
\\( P(\mathbf k) = \sum\_\nu f\_\nu c\_\nu c\_\nu^\dagger \\) を使って

$$
\rho(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\\, e^{-i\mathbf k\cdot\mathbf L}\\, P(\mathbf k)
$$

とも書けます。これは Bloch 和 \\( P(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L}\rho(\mathbf L) \\) の逆変換そのもの（k 点の和が厳密な場合）です。
したがって、実空間の行列の積は k 空間の積に対応し、\\( \rho \\)、\\( S \\)、\\( H \\) の間の関係（下の和則と冪等性）を実空間の畳み込みで検査できます。

### 非共線の場合〔規約〕

非共線計算では、スピノル固有ベクトル（前半が上向き、後半が下向きスピンの係数）から

$$
\rho^{\sigma\sigma'}\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\nu} f\_{\nu}(\mathbf k)\\,
e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{\sigma I\alpha,\nu}(\mathbf k)}\\, c\_{\sigma' J\beta,\nu}(\mathbf k)
$$

を作ります。複素数で、既存実装は次の成分を保持します。

- \\( \mathrm{Re}\\,\rho^{\uparrow\uparrow} \\)、\\( \mathrm{Re}\\,\rho^{\downarrow\downarrow} \\)、\\( \mathrm{Re}\\,\rho^{\uparrow\downarrow} \\)、\\( \mathrm{Im}\\,\rho^{\uparrow\downarrow} \\) の 4 つ（常に）。
- \\( \mathrm{Im}\\,\rho^{\uparrow\uparrow} \\)、\\( \mathrm{Im}\\,\rho^{\downarrow\downarrow} \\)（スピン軌道相互作用、LDA+U、スピン拘束、Zeeman 項のいずれかがあるとき）。

\\( \rho^{\downarrow\uparrow} \\) は Hermite 性 \\( \rho^{\sigma'\sigma}\_{J\beta,I\alpha}(-\mathbf L) = \overline{\rho^{\sigma\sigma'}\_{I\alpha,J\beta}(\mathbf L)} \\) から得られます。
これは [k 点と Bloch 和](kpoints-bloch.md) が「自然な形」として挙げた式と同じで、既存実装もこの形です。
対角スピン成分の虚部は、時間反転が破れている場合だけ 0 でなくなり、格子上の密度には寄与しません（実の基底関数の積は \\( a\leftrightarrow b \\) について対称）が、
軌道磁気モーメントや LDA+U の占有行列に使われます。

### エネルギー重み付き密度行列〔規約〕

力の計算（重なり行列の微分の項）のため、同じ形で固有値を掛けた

$$
E^{\sigma}\_{I\alpha, J\beta}(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k} \sum\_{\nu} f\_{\nu\sigma}(\mathbf k)\\, \varepsilon\_{\nu\sigma}(\mathbf k)\\,
\mathrm{Re}\left[ e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\nu\sigma}(\mathbf k)}\\, c\_{J\beta,\nu\sigma}(\mathbf k) \right]
$$

を同時に作ります。使い方は [全エネルギーと力](energy-force.md) にあります。\\( E \\) はミキシングしません（常に最後の固有値問題の出力です）。

### 和則〔設計〕

定義から直ちに次が成り立ちます（\\( \sum\_{\mathbf L}\sum\_{I\alpha J\beta} \\) を \\( \sum \\) と略記）。テストの不変量に使います。

$$
g\_s \sum\_{\sigma} \sum \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, S\_{I\alpha,J\beta}(\mathbf L) = N\_{\mathrm e},\qquad
g\_s \sum\_{\sigma} \sum \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, H^{\sigma}\_{I\alpha,J\beta}(\mathbf L) = U\_{\mathrm{ele}},\qquad
g\_s \sum\_{\sigma} \sum E^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\, S\_{I\alpha,J\beta}(\mathbf L) = U\_{\mathrm{ele}}
$$

1 番目は \\( c^\dagger S(\mathbf k) c = 1 \\)、2 番目は \\( c^\dagger H(\mathbf k) c = \varepsilon \\) から従います。
さらに \\( T = 0 \\)（占有数が 0 か 1）の絶縁体では、k 空間の \\( P S P = P \\) から、実空間で

$$
\sum\_{\mathbf L\_1 + \mathbf L\_2 + \mathbf L\_3 = \mathbf L} \rho(\mathbf L\_1)\\, S(\mathbf L\_2)\\, \rho(\mathbf L\_3) = \rho(\mathbf L)
$$

が成り立ちます。k 点が \\( N\_1\times N\_2\times N\_3 \\) の Γ 中心メッシュなら、\\( \mathbf L \\) を \\( N\_i \\) 周期で同一視した巡回畳み込みとして厳密に成り立ちます
（中心をずらしたメッシュでは \\( e^{i\mathbf k\cdot N\_i\mathbf a\_i} \ne 1 \\) なので巡回にならない点に注意）。
ただし \\( \rho \\) は疎ブロック形式の範囲（\\( \lVert\boldsymbol\tau\rVert < r\_{c,I} + r\_{c,J} \\)）でしか作らないので、
実際の計算の \\( \rho \\) では、範囲外の成分を落とした分だけ等式が破れます。

## 3. 格子上の電子密度〔公開〕〔規約〕

〔公開〕技術ノート 1（Eq. 4 前後）のとおり、

$$
n^{\sigma}(\mathbf r) = \sum\_{I}\sum\_{(J,\mathbf L)}\sum\_{\alpha\beta} \rho^{\sigma}\_{I\alpha,J\beta}(\mathbf L)\\,
\phi\_{I\alpha}(\mathbf r - \mathbf R\_I)\\, \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L)
$$

を実空間格子の各点で評価します（基底関数は実数）。格子、周期像の扱い、基底関数の格子上の値の持ち方は
[実空間グリッド・Poisson・交換相関](grid-poisson-xc.md) を参照してください。

〔規約〕

- 非スピン偏極では \\( n^{\uparrow} = n^{\downarrow} \\) で、それぞれ 1 スピン分の密度行列から作ります（全密度は和）。
- 非共線では、\\( n^{\uparrow\uparrow} \\)、\\( n^{\downarrow\downarrow} \\)、\\( \mathrm{Re}\\,n^{\uparrow\downarrow} \\)、\\( \mathrm{Im}\\,n^{\uparrow\downarrow} \\) の 4 つを格子上に持ちます。
  交換相関の計算では各点で 2×2 の密度行列を対角化して局所的なスピン軸に回します（詳細は格子の章）。
- 原子密度の重ね合わせ \\( n\_{\mathrm{atom}}(\mathbf r) \\)（初期密度と Hartree ポテンシャルの基準）は、スピンごとに半分ずつの形でも別に持ちます。

## 4. 電荷ミキシング

### 記法〔設計〕

反復 \\( m \\) の入力を \\( x^{\mathrm{in}}\_m \\)、それから得た出力を \\( x^{\mathrm{out}}\_m \\)、残差を \\( R\_m = x^{\mathrm{out}}\_m - x^{\mathrm{in}}\_m \\) と書きます。
\\( x \\) は方式によって、密度行列 \\( \rho(\mathbf L) \\)、格子上の密度の Fourier 成分 \\( n(\mathbf G) \\)、ポテンシャルの Fourier 成分 \\( V(\mathbf G) \\)、
または Hamiltonian \\( H(\mathbf L) \\) です。履歴は新しい順に \\( m = n, n-1, \dots \\) と並べます。

格子上の量の Fourier 成分は規格化なしの離散 Fourier 変換

$$
\tilde n(\mathbf G) = \sum\_{\mathbf p} n(\mathbf r\_{\mathbf p})\\, e^{-i\mathbf G\cdot\mathbf r\_{\mathbf p}},\qquad
n(\mathbf r\_{\mathbf p}) = \frac{1}{N\_{\mathrm g}} \sum\_{\mathbf G} \tilde n(\mathbf G)\\, e^{+i\mathbf G\cdot\mathbf r\_{\mathbf p}}
$$

で定めます（\\( N\_{\mathrm g} = N\_1 N\_2 N\_3 \\) は格子点の総数、\\( \mathbf G \\) は \\( -N\_i/2 \\) 付近から \\( N\_i/2 \\) 付近までの整数の組で表した逆格子ベクトル）。
〔公開〕技術ノート 12 の Eq. 7 は前向き変換に \\( 1/N\_{\mathrm g} \\) を付けていますが、ミキシングは線形なので結果は変わりません。
残差ノルムの値だけが規格化によって変わるので、下の〔規約〕の正規化はこの定義に対するものです。

### 方式の一覧〔公開〕〔規約〕

〔公開〕マニュアル 16.1 節と 8.2 節の `scf.Mixing.Type` の説明によります。〔規約〕の列は既存実装で確認した事実です。

| キーワード値 | 混ぜる量 | 空間 | 備考 |
|---|---|---|---|
| `Simple` | 密度行列 | 実空間（行列要素） | 〔規約〕既定値 |
| `RMM-DIIS` | 密度行列 | 実空間 | 〔規約〕最初の数反復は Hamiltonian を混ぜる |
| `GR-Pulay` | 密度行列 | 実空間 | 〔規約〕Pulay 段階では混合と無混合を交互に行う |
| `Kerker` | スピンごとの密度 \\( \tilde n^{\sigma}(\mathbf G) \\) | 逆空間 | |
| `RMM-DIISK` | スピンごとの密度 \\( \tilde n^{\sigma}(\mathbf G) \\) | 逆空間 | マニュアルが第一候補として推奨 |
| `RMM-DIISV` | 密度に依存する局所ポテンシャル \\( \tilde V^{\sigma}(\mathbf G) \\) | 逆空間 | |
| `RMM-DIISH` | Hamiltonian の行列要素 | 実空間 | LDA+U・拘束計算向け |
| `RMM-ADIIS` | 密度行列 | 実空間 | 〔規約〕実装にはあるがマニュアルに記載なし（〔未確定〕） |

〔規約〕

- `RMM-DIISH` と `RMM-DIISV` は軌道最適化（`orbitalOpt.*`）と併用できません。`RMM-DIIS` では `scf.Mixing.StartPulay` が 4 以上でなければなりません。
- 逆空間の方式で混ぜるのは、全密度そのもの（スピンごと）です。技術ノート 12 は原子密度の和からの差 \\( \delta n \\) で書いていますが、
  原子密度は SCF 中に変わらず、ミキシングの係数の和は常に 1 なので、どちらを混ぜても同じ結果になります。
- 非共線の逆空間の方式は、\\( \tilde n^{\uparrow\uparrow} \\)、\\( \tilde n^{\downarrow\downarrow} \\)、\\( \tilde n^{\uparrow\downarrow} \\)（複素）の 3 成分を、同じ係数で混ぜます。
- `RMM-DIISV` の \\( V^{\sigma} \\) は、Hartree ポテンシャル（原子密度の和からの差による部分）と交換相関ポテンシャルの和です
  （中性原子ポテンシャルのように密度に依存しない部分は含みません）。
- 密度行列を混ぜる方式では、混ぜた \\( \rho^{\mathrm{in}} \\) から格子上の \\( n^{\mathrm{in}} \\) を作ります。逆空間の方式では、
  出力の \\( \rho^{\mathrm{out}} \\) から \\( n^{\mathrm{out}} \\) を作って Fourier 変換し、混ぜてから逆変換して \\( n^{\mathrm{in}} \\) にします（\\( \rho \\) 自体は混ぜません）。

### 単純混合〔公開〕〔規約〕

〔公開〕技術ノート 12 の Eq. 1 は、入力と出力の線形結合で次の入力を作る単純混合です。

〔規約〕既存実装の重み \\( w \\)（`scf.*.Mixing.Weight`）は **出力側の割合** です：

$$
x^{\mathrm{in}}\_{n+1} = w\\, x^{\mathrm{out}}\_n + (1 - w)\\, x^{\mathrm{in}}\_n = x^{\mathrm{in}}\_n + w R\_n
$$

技術ノート 12 の Eq. 1 の \\( \alpha \\) は入力側の割合なので、\\( \alpha = 1 - w \\) と読み替えます。
密度行列の単純混合は、保持しているすべての成分（スピン成分、非共線の実部・虚部）に同じ \\( w \\) を使います。

- 最初の MD ステップの最初の反復では \\( w = 1 \\)（出力をそのまま次の入力にする）。
- MD・構造最適化の 2 ステップ目以降は、各 SCF の最初の重みを `scf.Max.Mixing.Weight` から始めます。

### 残差ノルム〔公開〕〔規約〕

〔公開〕マニュアル 16.1 節によると、各反復の標準出力の `NormRD` が残差の大きさです。

〔規約〕既存実装は残差の 2 乗ノルム \\( \mathcal N \\) を保持し、表示する NormRD はその平方根 \\( \sqrt{\mathcal N} \\) です。
\\( \mathcal N \\) の定義は方式ごとに異なります。

| 方式・段階 | \\( \mathcal N \\) |
|---|---|
| 密度行列（`Simple`、`RMM-DIIS`、`GR-Pulay`） | \\( \sum\_{\text{成分}}\sum\_{(I,J,\mathbf L)}\sum\_{\alpha\beta} \left(R\_{I\alpha,J\beta}(\mathbf L)\right)^2 \\)（正規化なし） |
| `Kerker`、および RMM-DIISK の Pulay 開始前 | \\( \frac{1}{N\_{\mathrm g}} \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.2}(\mathbf G)\\, \lvert\tilde R^{\sigma}(\mathbf G)\rvert^2 \\) |
| RMM-DIISK の Pulay 段階 | \\( \frac{1}{N\_{\mathrm g} N\_s} \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.1}(\mathbf G)\\, \lvert\tilde R^{\sigma}(\mathbf G)\rvert^2 \\) |
| RMM-DIISV | \\( \frac{1}{10^{5} N\_{\mathrm g} N\_s} \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.1}(\mathbf G)\\, \lvert\tilde R\_V^{\sigma}(\mathbf G)\rvert^2 \\) |
| Hamiltonian の単純混合 | \\( \frac{1}{N\_{\mathrm{atom}}} \sum (R\_H)^2 \\) |
| RMM-DIISH | \\( \frac{1}{N\_{\mathrm{atom}}} \sum m\_{I\alpha} (R\_H)^2 \\)（\\( m \\) は下の計量） |

ここで \\( N\_s \\) は逆空間で混ぜるスピン成分の数（1、2、3）、\\( N\_{\mathrm{atom}} \\) は単位胞の原子数、
\\( K\_{\gamma} \\) は Kerker 型の計量（次節）です。\\( G = 0 \\) 成分の計量の値は、Kerker と RMM-DIISV で 1、RMM-DIISK の Pulay 段階で 0.5 です。
各 SCF の最初の反復では \\( \mathcal N \\) を 1 とみなします（初期値は 100）。

\\( \mathcal N \\) は方式によって単位も正規化も異なるので、方式をまたいで NormRD の値を比べることはできません。

### 重みの自動調整〔規約〕

単純混合（密度行列、Kerker、Hamiltonian、および RMM-DIISK・RMM-DIISV で Pulay を行わない反復）の重み \\( w \\) は、
残差ノルムの変化に応じて毎反復更新します。現在と 1 つ前の 2 乗ノルムを \\( \mathcal N\_0 \\)、\\( \mathcal N\_1 \\) として、

$$
w \leftarrow
\begin{cases}
w\\, \dfrac{\mathcal N\_1}{\max(\mathcal N\_1 - \mathcal N\_0,\\, 10^{-9})} & (\mathcal N\_0 < \mathcal N\_1\text{：残差が減った})\\\\[2mm]
w\\, \dfrac{\mathcal N\_1}{\mathcal N\_1 + \mathcal N\_0} & (\mathcal N\_0 > \mathcal N\_1\text{：残差が増えた})\\\\[2mm]
w & (\mathcal N\_0 = \mathcal N\_1)
\end{cases}
$$

とし、\\( [w\_{\min}, w\_{\max}] \\) に切り詰めます（`scf.Min.Mixing.Weight`、`scf.Max.Mixing.Weight`）。

- 残差の減りが遅いほど（\\( \mathcal N\_0 \approx \mathcal N\_1 \\)）重みを大きくし、残差が増えたら少なくとも半分にします。
- 初期値は `scf.Init.Mixing.Weight` です。
- 上限に張り付いたまま残差が減る反復が非常に多く（1000 回）続いた場合だけ、上限を 2 倍（最大 0.7）にする仕組みもありますが、実用上は働きません。

RMM-DIISV の Pulay 開始前の段階（ポテンシャルの Kerker 混合）だけは別の規則で、直近 3 回のバンドエネルギー \\( U\_0, U\_1, U\_2 \\) の変化の向きと残差の増減で

| \\( U \\) の変化の向き | 残差 | 新しい \\( w \\) |
|---|---|---|
| 同じ（\\( \mathrm{sgn}(U\_0 - U\_1) = \mathrm{sgn}(U\_1 - U\_2) \\)） | 減少 | \\( 1.6\\, w \\) |
| 同じ | 増加 | \\( w / 1.6 \\) |
| 逆 | 減少 | \\( 1.2\\, w \\) |
| 逆 | 増加 | \\( w / 2 \\) |

とします（同じ範囲に切り詰め）。

### Pulay 型の係数（RMM-DIIS）〔公開〕〔規約〕

〔公開〕技術ノート 12 の Eq. 2–6（Pulay, CPL **73**, 393 (1980)）：直近 \\( p \\) 個の残差の線形結合
\\( \bar R = \sum\_m \alpha\_m R\_m \\) のノルムを、\\( \sum\_m \alpha\_m = 1 \\) の条件で最小にします。
Lagrange 乗数 \\( \lambda \\) を使うと、

$$
\begin{pmatrix} A & \mathbf 1 \\\\ \mathbf 1^{\mathsf T} & 0 \end{pmatrix}
\begin{pmatrix} \boldsymbol\alpha \\\\ \lambda/2 \end{pmatrix} = \begin{pmatrix} \mathbf 0 \\\\ 1 \end{pmatrix},\qquad
A\_{mm'} = \langle R\_m \vert R\_{m'} \rangle
$$

で、閉じた形では \\( \boldsymbol\alpha = A^{-1}\mathbf 1 / (\mathbf 1^{\mathsf T}A^{-1}\mathbf 1) \\) です。

〔規約〕既存実装は、得られた \\( \alpha \\) で入力の線形結合に残差の補正を加えます：

$$
x^{\mathrm{in}}\_{n+1} = \sum\_{m} \alpha\_m x^{\mathrm{in}}\_m + \beta\\, \bar R
$$

技術ノート 12 の Eq. 6 は \\( \beta = 0 \\) の形ですが、既存実装は \\( \beta \\) を現在の 2 乗ノルム \\( \mathcal N \\) に応じて次の表から選びます。

| 方式 | \\( \mathcal N \\) による \\( \beta \\) |
|---|---|
| RMM-DIIS（密度行列）、RMM-DIISH | \\( \mathcal N \ge 10^{-1} \\)：0.5、\\( \ge 10^{-2} \\)：0.6、\\( \ge 10^{-3} \\)：0.7、\\( \ge 10^{-4} \\)：0.8、それ未満：1.0 |
| RMM-DIISK、RMM-DIISV | \\( \mathcal N \ge 10^{-2} \\)：0.5、\\( \ge 10^{-4} \\)：0.6、\\( \ge 10^{-6} \\)：0.7、\\( \ge 10^{-8} \\)：0.8、\\( \ge 10^{-10} \\)：0.9、それ未満：1.0 |
| GR-Pulay | \\( \mathcal N \ge 10^{-2} \\)：0.2、\\( \ge 10^{-4} \\)：0.3、\\( \ge 10^{-6} \\)：0.4、\\( \ge 10^{-8} \\)：0.5、それ未満：1.0 |

- 内積 \\( \langle R\_m\vert R\_{m'}\rangle \\) は、密度行列ではすべての成分の積の和、逆空間と Hamiltonian では下の計量つきの和です。
- 1 次方程式は小さい（履歴の数 + 1 次元）ので、直接解きます。解に NaN・無限大が出た場合は、
  最新の 2 つだけを使う \\( (\alpha\_n, \alpha\_{n-1}) = (0.1, 0.9) \\)（RMM-DIISH は \\( (0.05, 0.95) \\)）に置き換えます。
- 履歴の長さは `scf.Mixing.History`（既定 5）が上限です。Pulay 段階に入ると 2 本から始め、1 反復ごとに 1 本ずつ増やします。

〔設計〕\\( A \\) は残差が一次従属に近づくと特異に近づきます。再実装では、対角成分で規格化してから解く、
最小固有値が閾値を下回ったら古い履歴を捨てる、などの安定化を入れてかまいません（既存実装の NaN 時の置き換えと同じ目的）。

### RMM-DIIS（密度行列）の段階〔規約〕

`scf.Mixing.StartPulay` を \\( P \\)（既定 6）、SCF の反復番号を \\( n \\) とします。

1. \\( n \le \lfloor P/2\rfloor \\)：密度行列は混ぜず（出力をそのまま保持）、代わりに Hamiltonian の単純混合を、固有値問題を解く前に行います：
   \\( H^{\mathrm{in}}\_{n} = w\\, H[n^{\mathrm{in}}\_n] + (1-w)\\, H^{\mathrm{in}}\_{n-1} \\)。重みは自動調整します（残差は \\( H[n^{\mathrm{in}}\_n] - H^{\mathrm{in}}\_{n-1} \\)）。
   restart ファイルから始めた場合と軌道最適化を使う場合は、この段階を省きます。
2. \\( \lfloor P/2\rfloor < n < P \\)：密度行列の単純混合（重みは自動調整）。
3. \\( n = P \\)：重みを固定した（調整しない）単純混合。残差の履歴を Pulay 用に並べ直します。
4. \\( n > P \\)：Pulay 混合（履歴 2 本から）。

GR-Pulay も同じく \\( n < P \\) で単純混合、\\( n = P \\) で重み固定の単純混合を行い、\\( n > P \\) から上の交互の手順に入ります。

### GR-Pulay〔規約〕

Bowler と Gillan の「保証付き残差減少」型の Pulay 法（Chem. Phys. Lett. **325**, 473 (2000)）に相当する方式です。
Pulay 開始前は密度行列の単純混合です（RMM-DIIS のような Hamiltonian の混合の段階はありません。最初の反復は \\( w = 1 \\)）。Pulay 段階では、反復を 2 つ 1 組にして交互に次を行います。

- 奇数番目：現在の出力 \\( \rho^{\mathrm{out}} \\) を保存し、その残差を履歴に入れ、**混ぜずに** 出力をそのまま次の入力にします。
- 偶数番目：入力は前の出力 \\( \tilde\rho = \rho^{\mathrm{out}}\_{\text{前}} \\) なので、残差 \\( R\_0 = \rho^{\mathrm{out}}[\tilde\rho] - \tilde\rho \\) を作り、
  \\( \\{(\tilde\rho, R\_0)\\} \\) と履歴 \\( \\{(x^{\mathrm{in}}\_m, R\_m)\\} \\) で Pulay 係数を求めて、\\( x^{\mathrm{in}} = \alpha\_0 \tilde\rho + \sum\_m \alpha\_m x^{\mathrm{in}}\_m + \beta\bar R \\) とします。
  係数は \\( A \\) を \\( A\_{00} \\) で割ってから閉じた形で求めます。

既存実装の GR-Pulay は非共線の虚部成分（スピン軌道相互作用などがあるときの対角スピン成分の虚部）を混ぜません。

### Kerker 混合〔公開〕〔規約〕

〔公開〕技術ノート 12 の Eq. 7–10（Kerker, PRB **23**, 3082 (1981)）：逆空間で成分ごとに重みを変えた単純混合で、
重み \\( g(\mathbf G) = \lvert\mathbf G\rvert^2/(\lvert\mathbf G\rvert^2 + q\_0^2) \\) は長波長成分（電荷のスロッシングの主因）の更新を抑えます。
\\( q\_0 = \gamma\lvert\mathbf G\_{\min}\rvert \\) で、\\( \gamma \\) が `scf.Kerker.factor` です。

〔規約〕既存実装の重みは、長波長で 0 にならない形です。最短の逆格子ベクトルの長さを \\( b\_{\min} = \min\_i\lVert\mathbf b\_i\rVert \\)
（\\( \mathbf b\_i\cdot\mathbf a\_j = 2\pi\delta\_{ij} \\)）、\\( q\_0 = \gamma\\, b\_{\min} \\) として、計量

$$
K\_{\eta}(\mathbf G) = \frac{G^2 + q\_0^2}{G^2 + (\eta q\_0)^2}\quad (\mathbf G \ne 0)
$$

を定義すると、Kerker 混合は

$$
\tilde n^{\mathrm{in}}\_{n+1}(\mathbf G) = \tilde n^{\mathrm{in}}\_n(\mathbf G) + \frac{w}{\sqrt{K\_{0.2}(\mathbf G)}}\\, \tilde R\_n(\mathbf G),\qquad K\_{0.2}(0) = 1
$$

です。混合の割合は大きな \\( G \\) で \\( w \\)、\\( G \to 0 \\) で \\( 0.2\\, w \\)、\\( G = 0 \\) で \\( w \\) です。

- 重み \\( w \\) は前述の規則で自動調整します。
- 最初の MD ステップの最初の反復では、既存実装は **入力をそのまま保ちます**（出力側の割合 0。この反復の \\( K \\) は 1）。
  結果として 2 回目の反復は 1 回目と同じ Hamiltonian になります。〔設計〕再実装では、この反復を \\( w = 1 \\) の単純混合
  （出力をそのまま採用）にするか、Kerker 混合にしてかまいません。SCF の経過が 1 反復ずれるだけです。
- 長波長で重みが \\( 0.2w \\) に下がるだけなので、大きな金属の単位胞では、教科書の Kerker 混合（\\( G \to 0 \\) で重み 0）ほど
  スロッシングを抑えません（`tools/check_density_mixing.py`：Thomas–Fermi 型の応答 \\( \chi \propto q\_{\mathrm{TF}}^2/G^2 \\) では、
  最長波長の成分の増幅率は \\( 1 - 0.2w(1 + q\_{\mathrm{TF}}^2/G\_{\min}^2) \\) で、\\( G\_{\min} \ll q\_{\mathrm{TF}} \\) なら \\( w \\) を十分小さくしないと発散）。
  実際には重みの自動調整が \\( w \\) を下げて対応します。

### RMM-DIISK〔公開〕〔規約〕

〔公開〕技術ノート 12 の Eq. 11–12：逆空間の密度に Pulay 法を使い、内積に Kerker の計量 \\( 1/g(\mathbf G) \\) を使います。
\\( 1/g = (G^2 + q\_0^2)/G^2 \\) は長波長の残差を重く見ます。

〔規約〕

- 計量は \\( G \to 0 \\) で発散しないように正則化した \\( K\_{0.1}(\mathbf G) \\)、\\( G = 0 \\) では 0.5 です：
  \\( \langle R\_m\vert R\_{m'}\rangle = \sum\_{\sigma}\sum\_{\mathbf G} K\_{0.1}(\mathbf G)\\, \mathrm{Re}\left[\overline{\tilde R\_m^{\sigma}(\mathbf G)}\\, \tilde R\_{m'}^{\sigma}(\mathbf G)\right] \\)。
- Pulay 開始前（\\( n < P \\)）は上の Kerker 混合です。この間の残差も履歴に入るので、Pulay 段階の最初から履歴が使えます。
- Pulay 段階の補正項は計量で割り戻します：
  \\( \tilde n^{\mathrm{in}} = \sum\_m \alpha\_m \tilde n^{\mathrm{in}}\_m + \beta\\, \bar R' \\) で、
  \\( \mathcal N \le 1 \\) なら \\( \bar R' = \sum\_m\alpha\_m \tilde R\_m \\)、\\( \mathcal N > 1 \\) なら \\( \bar R' = \sum\_m\alpha\_m \tilde R\_m / \sqrt{K\_{0.1}} \\)（長波長をさらに抑える）。
- `scf.Mixing.EveryPulay` を \\( E \\)（既定 1）とすると、Pulay 段階では \\( E \\) 回に 1 回だけ Pulay 混合を行い、その間は
  \\( \tilde n^{\mathrm{in}}\_{n+1} = \tilde n^{\mathrm{in}}\_n + w\\, \tilde R\_n / K\_{0.1}(\mathbf G) \\)（\\( G = 0 \\) では \\( 2w \\)）の Kerker 型混合を行います。
  Pulay 混合を行うたびに、次の Kerker 型混合のために \\( w \\) を \\( \max(0.1\\, w\_{\max},\\, w\_{\min}) \\) に戻します。
- 履歴の数え方には癖があります：Pulay 段階の最初の反復では \\( \min(P-1, h) \\) 本（\\( P \le 6 \\) のとき）または 2 本（\\( P > 6 \\) のとき）を使い、
  その後は 1 反復ごとに 1 本ずつ増やして `scf.Mixing.History`（\\( h \\)）を上限とします。
  （\\( P \le 6 \\) では 2 反復目に本数が減るので、再実装では \\( P \\) によらず 2 本から始めて増やす形で十分です。）

〔公開〕マニュアル 16.1 節は `scf.Mixing.EveryPulay` が RMM-DIISK・RMM-DIISV・RMM-DIISH に効くとしていますが、
既存実装で効くのは RMM-DIISK と RMM-DIISV だけです（RMM-DIISH は Pulay 段階で毎回 Pulay 混合）。

### RMM-DIISV〔規約〕

RMM-DIISK と同じ手順を、密度の代わりにポテンシャル \\( \tilde V^{\sigma}(\mathbf G) \\) に適用します。Hamiltonian を組み立てる直前に混ぜます。

- 計量は \\( K\_{0.1} \\)、ただし \\( G = 0 \\) では 1。2 乗ノルムの正規化は前の表のとおり（\\( 10^5 \\) で割る）。
- Pulay 開始前の Kerker 型混合は \\( \tilde V^{\mathrm{in}}\_{n+1} = \tilde V^{\mathrm{in}}\_n + (w/\sqrt{K\_{0.1}})\\, \tilde R\_V \\) で、重みはバンドエネルギーの傾向を使う別の規則（前述）で調整します。
- 最初の反復では新しいポテンシャルをそのまま使います。
- Pulay 段階の \\( \beta \\) の表、`scf.Mixing.EveryPulay`、補正項の割り戻しは RMM-DIISK と同じです。

### RMM-DIISH〔規約〕

Hamiltonian の行列要素 \\( H(\mathbf L) \\) を、固有値問題を解く直前に混ぜます。
LDA+U や拘束の項は密度行列（占有行列）を通じて Hamiltonian に入るので、密度だけを混ぜるより安定になります（マニュアル 16.1 節の推奨）。

- 入力は、前の反復で混ぜた Hamiltonian \\( H^{\mathrm{in}}\_{n-1} \\)、出力は今回の入力の密度から組み立てた \\( H[n^{\mathrm{in}}\_n] \\) で、残差はその差です。
- Pulay 開始前（\\( n < P \\)）は Hamiltonian の単純混合（重みは自動調整）。最初の反復は \\( w = 1 \\)。
- Pulay 段階の内積は、原子 \\( I \\) の軌道 \\( \alpha \\)（行の添字）ごとの計量

$$
m\_{I\alpha} = \frac{5}{d^2 + 5},\qquad d = \left\lvert H^{\uparrow\uparrow}\_{I\alpha,I\alpha}(\mathbf 0) - \mu \right\rvert\ \ (\text{Hartree})
$$

  を掛けた和です（\\( H \\) は直近の入力、\\( \mu \\) は化学ポテンシャル）。Fermi 準位から遠いエネルギーの軌道ほど軽く見ます。
- 非共線では \\( H^{\uparrow\uparrow} \\)、\\( H^{\downarrow\downarrow} \\)、\\( \mathrm{Re}\\,H^{\uparrow\downarrow} \\)、\\( \mathrm{Im}\\,H^{\uparrow\downarrow} \\) と、
  スピン軌道相互作用などによる虚部の行列のすべてを同じ係数で混ぜます。

### Kerker 因子の自動決定〔公開〕〔規約〕

〔公開〕マニュアル 16.2 節：`scf.Kerker.factor` を指定しない場合、逆格子ベクトルの長さの平均と非等方性から自動で決めます（同節 Eq. 21–23）。
マニュアルの式は逆格子ベクトルの長さの 2 乗で平均と差を作り、係数 4 を掛ける形に見えますが、表示上の分数の構造が読み取りにくく、下の実装と細部が一致しません。

〔規約〕既存実装（Kerker、RMM-DIISK、RMM-DIISV で、指定がない場合）は、逆格子ベクトルの長さ \\( b\_i = \lVert\mathbf b\_i\rVert \\) から

$$
\bar b = \frac{b\_1 + b\_2 + b\_3}{3},\qquad
\Delta b = \frac{\lvert b\_2 - b\_1\rvert + \lvert b\_3 - b\_2\rvert + \lvert b\_1 - b\_3\rvert}{3},\qquad
\gamma = \frac{0.5}{b\_{\min}}\left(\frac{2\Delta b}{\bar b} + 1\right)
$$

とします。したがって \\( q\_0 = \gamma b\_{\min} = 0.5\\,(2\Delta b/\bar b + 1) \\) bohr⁻¹ で、単位胞の大きさによらず、
等方的な単位胞で 0.5 bohr⁻¹、1 方向だけ極端に長い単位胞（細長い単位胞、逆格子ベクトルが 1 本だけ短い）で 1.5 bohr⁻¹ に近づき、
1 方向だけ極端に短い単位胞（逆格子ベクトルが 1 本だけ長い）で最大 2.5 bohr⁻¹ です
（例：5×5×40 bohr の直方体で 1.32 bohr⁻¹）。
指定した場合は \\( q\_0 = \gamma\\, b\_{\min} \\) で、\\( b\_{\min} \\) の分だけ単位胞の大きさに依存します。

### パラメータと既定値〔公開〕〔規約〕

〔公開〕の既定値はマニュアル 8.2 節の記載、〔規約〕は既存実装で確認した値です（両者は一致しています）。

| キーワード | 既定値 | 意味 |
|---|---|---|
| `scf.Mixing.Type` | `Simple`〔規約〕 | 方式 |
| `scf.Init.Mixing.Weight` | 0.3〔公開〕 | 単純混合の重みの初期値（出力側の割合） |
| `scf.Min.Mixing.Weight` | 0.001〔公開〕 | 重みの下限 |
| `scf.Max.Mixing.Weight` | 0.4〔公開〕 | 重みの上限。MD・最適化の 2 ステップ目以降の初期値 |
| `scf.Kerker.factor` | 自動〔公開〕 | \\( \gamma \\)。上の自動決定 |
| `scf.Mixing.History` | 5〔規約〕 | Pulay 型の履歴の上限。マニュアルは 30 前後を推奨 |
| `scf.Mixing.StartPulay` | 6〔公開〕 | Pulay 混合を始める反復 |
| `scf.Mixing.EveryPulay` | 1〔公開〕 | RMM-DIISK・RMM-DIISV で Pulay を行う間隔 |
| `scf.Mixing.Control.Temp` | off〔規約〕 | 電子温度の自動制御 |
| `scf.criterion` | \\( 10^{-6} \\) Hartree〔公開〕 | 収束判定の閾値 |
| `scf.maxIter` | 40〔公開〕 | 反復の上限 |
| `scf.ElectronicTemperature` | 300 K〔公開〕 | 電子温度 |
| `scf.system.charge` | 0〔公開〕 | 余分な電荷 |
| `scf.ExtCharge.History` | 3〔規約〕 | 密度の外挿に使う履歴の数 |

### SCF の途中での調整〔公開〕〔規約〕

〔公開〕マニュアル 16.3 節：計算の途中で、作業ディレクトリに `System.Name` に `_SCF_keywords` を付けた名前のファイルを置くと、
SCF の各反復の終わりにそのファイルを読み、`scf.maxIter`、`scf.Min.Mixing.Weight`、`scf.Max.Mixing.Weight`、`scf.Kerker.factor`、
`scf.Mixing.StartPulay` を更新します。`scf.maxIter` を負にすると計算を止めます。

〔規約〕既存実装はこのファイルから `scf.criterion` も読みます。ファイルは読んだ後も残り、毎反復読み直されます。

〔規約〕MD・構造最適化の 2 ステップ目以降で、前のステップの密度から再開できた場合、
`scf.Mixing.StartPulay` より前の反復で \\( \mathcal N < 0.01 \\) になったら、次の反復から Pulay 混合を始めます（そのステップの SCF の間だけ）。

## 5. 収束判定〔公開〕〔規約〕

〔公開〕マニュアル 8.2 節と 16.1 節：`scf.criterion`（Hartree 単位、既定 \\( 10^{-6} \\)）は、反復間のバンドエネルギーの変化
\\( \mathrm{d}U\_{\mathrm{ele}} \\) に対する閾値です。`scf.maxIter`（既定 40）回で、収束していなくても SCF を打ち切ります。

〔規約〕既存実装の判定は、バンドエネルギーの変化だけでなく残差ノルムも見ます。閾値を \\( c \\) として、

$$
n \ge 2,\qquad \left\lvert U\_{\mathrm{ele},n} - U\_{\mathrm{ele},n-1}\right\rvert < c,\qquad \mathcal N < 10\\, c
$$

をすべて満たしたら収束とします。\\( \mathcal N \\) は直近の 2 乗ノルム（NormRD の 2 乗、方式ごとの正規化は前の表）で、
エネルギーの次元をもつ \\( c \\) と次元の異なる量を比べている点に注意してください。

- 軌道最適化を使う場合は、残差ノルムの条件を使いません。
- 打ち切られた場合は非収束として扱い、MD の次のステップでは、\\( \mathcal N > 0.3 \\) なら前のステップの密度を再利用しません。

〔設計〕再実装では、同じ判定（エネルギーと残差ノルムの両方）を既定とし、残差ノルムの正規化は方式ごとに上の表に合わせます。
判定を変える場合は、全エネルギーの許容誤差（\\( 10^{-6} \\) Hartree 程度）が満たされることをテストで確認します。

## 6. 初期密度と再開〔公開〕〔規約〕

〔規約〕restart ファイルを使わない最初の SCF では、入力の密度は原子密度の重ね合わせです。

- 原子 \\( I \\) の密度は、擬原子の価電子密度（球対称、積分が \\( Z\_I^{\mathrm{val}} \\)）を、入力の原子座標の欄で与えた
  上向き・下向きの電子数 \\( N\_I^{\uparrow} \\)、\\( N\_I^{\downarrow} \\) を使って \\( (N\_I^{\uparrow} + N\_I^{\downarrow})/Z\_I^{\mathrm{val}} \\) 倍し、
  上向きに \\( N\_I^{\uparrow}/(N\_I^{\uparrow}+N\_I^{\downarrow}) \\)、下向きに残りを割り振ったものです（非偏極では半分ずつ）。
  したがって初期密度の電子数は入力の \\( \sum\_I (N\_I^{\uparrow}+N\_I^{\downarrow}) \\) で、`scf.system.charge` はこの段階では使いません
  （最初の固有値問題から正しい電子数になります）。
- Hartree ポテンシャルの基準に使う原子密度の和 \\( n\_{\mathrm{atom}} \\) は、この拡大縮小をしない中性原子の密度の和です。
- 非共線では、入力で与えた各原子のスピンの向き（角度）に沿って、2×2 の密度 \\( n^{\sigma\sigma'} \\) を作ります：
  原子ごとの全密度 \\( n \\) と磁化 \\( \mathbf m \\) から \\( n^{\uparrow\uparrow} = (n + m\_z)/2 \\)、\\( n^{\downarrow\downarrow} = (n - m\_z)/2 \\)、
  \\( \mathrm{Re}\\,n^{\uparrow\downarrow} = m\_x/2 \\)、\\( \mathrm{Im}\\,n^{\uparrow\downarrow} = -m\_y/2 \\)。
- 最初の反復は、この密度から Hamiltonian を作るので、密度行列の入力は必要ありません。

〔公開〕マニュアル 17.1 節：`scf.restart on` で、`System.Name_rst` ディレクトリに保存した前の計算の結果から SCF を再開します。

〔規約〕再開には、スピン成分ごとの格子上の密度の原子密度からの差 \\( \Delta n^{\sigma} = n^{\sigma} - n\_{\mathrm{atom}}/2 \\) と、Hamiltonian の行列要素を使います。

- 新しい入力の密度は \\( n^{\sigma} = n\_{\mathrm{atom}}/2 + \Delta n^{\sigma} \\)（\\( n\_{\mathrm{atom}} \\) は現在の原子位置で作り直したもの）です。
  格子の点数が保存時と異なる場合は再開しません。
- 再開した場合も、通常は保存した Hamiltonian を使わず、読み込んだ密度から Hamiltonian を作り直します
  （スピン軌道相互作用、LDA+U、拘束、Zeeman 項がある場合は、保存した Hamiltonian を最初の反復の Hamiltonian に使います）。
- MD・構造最適化では、2 ステップ目以降は自動的に前のステップの結果から再開します（`scf.restart alloff` で無効）。

## 7. MD・構造最適化での密度の外挿〔公開〕〔規約〕

〔公開〕技術ノート 8：原子座標を 1 本のベクトル \\( \mathbf x \\)（\\( 3N\_{\mathrm{atom}} \\) 次元）とみなし、
新しい座標を過去 \\( p \\) ステップの座標の線形結合 \\( \sum\_m \alpha\_m \mathbf x\_m \\) で最もよく近似する係数を求め（Eq. 1–3。\\( \sum\_m\alpha\_m = 1 \\) の条件つき最小 2 乗）、
同じ係数で密度を外挿します（Eq. 4）：

$$
n^{\mathrm{in}}\_{n+1} = n\_{\mathrm{atom}}[\mathbf x\_{n+1}] + \sum\_{m=n-p+1}^{n} \alpha\_m\\, \Delta n^{\mathrm{out}}\_m,\qquad
\Delta n^{\mathrm{out}}\_m = n^{\mathrm{out}}\_m - n\_{\mathrm{atom}}[\mathbf x\_m]
$$

技術ノートは \\( p = 3 \\) が多くの場合に最適としています。
係数は \\( \lVert \sum\_m\alpha\_m\mathbf x\_m - \mathbf x\_{n+1}\rVert^2 - \lambda(1 - \sum\_m\alpha\_m) \\) の停留条件

$$
\begin{pmatrix} X & \mathbf 1 \\\\ \mathbf 1^{\mathsf T} & 0 \end{pmatrix}
\begin{pmatrix} \boldsymbol\alpha \\\\ \lambda/2 \end{pmatrix} = \begin{pmatrix} \mathbf b \\\\ 1 \end{pmatrix},\qquad
X\_{mm'} = \mathbf x\_m\cdot\mathbf x\_{m'},\quad b\_m = \mathbf x\_m\cdot\mathbf x\_{n+1}
$$

から求まります（技術ノートの Eq. 3 は右辺の表示が崩れていますが、停留条件から右辺は \\( \mathbf b \\) です）。

〔規約〕既存実装の係数の決め方は次のとおりで、技術ノートと異なり **制約なし** の最小 2 乗です。

- 履歴の数は `scf.ExtCharge.History`（\\( p \\)、既定 3）。MD のステップ番号が \\( p \\) を超え、かつ \\( p > 1 \\) のときだけ外挿します。
  それまでは \\( (\alpha\_n, \alpha\_{n-1}, \dots) = (1, 0, \dots) \\)（前のステップの差の密度をそのまま使う）です。
- 正規方程式 \\( X\boldsymbol\alpha = \mathbf b \\) を、\\( X \\) の固有分解で解きます（固有値の絶対値が \\( 10^{-15} \\) 未満の方向は捨てる擬似逆行列）。
- 得られた係数のどれかの絶対値が 10 を超えるか、\\( \lvert\sum\_m\alpha\_m - 1\rvert > 10^{-3} \\) なら、外挿をやめて \\( (1, 0, \dots) \\) に戻します。
  （座標の絶対値が大きく、変位が小さいので、制約なしでも通常は和がほぼ 1 になります。）
- 外挿は格子上の密度（スピン成分ごと）に対して行い、Hamiltonian と密度行列は外挿しません。
- 原子が単位胞の外に出て分数座標を \\( [0, 1) \\) に戻した場合は、座標の履歴も同じだけずらして、外挿が連続になるようにします。

〔設計〕再実装では、技術ノートの制約つきの式を使い、既存実装と同じ安全装置（係数の大きさの上限と、失敗時の \\( (1, 0, \dots) \\)）を付けます。
等加速度運動の軌跡（座標が時刻の 2 次式）では、どちらの方法でも \\( \boldsymbol\alpha = (3, -3, 1) \\) が得られ、新しい座標が厳密に再現されます。

## 実装への要求〔設計〕

- 占有数・化学ポテンシャルは、固有値の配列（k 点、スピン、状態）と重み、電子数、温度だけを入力とする純粋な関数にします。
- 密度行列の構築は、疎ブロック形式の構造（原子対の集合）を受け取り、値だけを埋めます。\\( \rho \\) と \\( E \\) を同時に作ります。
- ミキシングは「混ぜる量のベクトル空間」（内積・計量・線形結合）を抽象化し、Simple／Pulay／Kerker の各方式を、
  密度行列・逆空間の密度・ポテンシャル・Hamiltonian に共通に適用できるようにします。
  履歴の長さ、\\( \beta \\) の表、重みの調整則、Pulay の開始反復は設定値として持ちます。
- 既存実装の SCF の経過（各反復の NormRD）を再現したい場合と、収束の速さを優先する場合とで、既定の規則を切り替えられるようにします。

## テスト

`tools/check_density_mixing.py` で、下の解析・不変量・線形模型の項目を確認しています。

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 2 準位系（\\( \varepsilon\_1 = -0.3 \\)、\\( \varepsilon\_2 = 0.1 \\) Hartree、非偏極）：\\( N\_{\mathrm e} = 3 \\) の \\( \mu \\) を、方程式 \\( 2(f\_1 + f\_2) = 3 \\) を別の求根法で解いた値と比較 | \\( \mu \\)：\\( 10^{-9} \\) Hartree |
| 解析 | 同じ系で \\( N\_{\mathrm e} = 2 \\)：電子数が合い、\\( \mu \\) がギャップ内にあること（中央との差は温度と停止条件しだい） | 電子数 \\( 10^{-11} \\) |
| 解析 | 二分法の停止条件 \\( \lvert N(\mu) - N\_{\mathrm e}\rvert < 10^{-12} \\) | 電子数 \\( 10^{-11} \\) |
| 不変量 | 電子数の和則 \\( g\_s\sum\rho S = N\_{\mathrm e} \\)（1 次元 2 軌道鎖、重なりあり、金属的温度と低温） | \\( 10^{-12} \\) |
| 不変量 | バンドエネルギーの和則 \\( g\_s\sum\rho H = g\_s\sum E S = U\_{\mathrm{ele}} \\) | \\( 10^{-12} \\) |
| 不変量 | Hermite 性 \\( \rho\_{J\beta,I\alpha}(-\mathbf L) = \rho\_{I\alpha,J\beta}(\mathbf L) \\)（Cluster と Band。O(N) 解法では成り立たない。非共線では複素共役とスピンの入れ替え） | 丸め誤差 |
| 不変量 | 絶縁体・低温での冪等性 \\( \rho S \rho = \rho \\)（Γ 中心メッシュ上の巡回畳み込み） | \\( 10^{-10} \\) |
| 不変量 | 非共線で \\( \rho^{\uparrow\downarrow} \\) を 0 にした系（共線の Hamiltonian を 2 つ並べたもの）の密度行列が、共線計算と一致 | \\( 10^{-12} \\) |
| 不変量 | 時間反転で \\( \pm\mathbf k \\) を重み 2 の 1 点にまとめても密度行列が変わらない | \\( 10^{-12} \\) |
| 線形模型 | \\( x^{\mathrm{out}} - x^\ast = J(x^{\mathrm{in}} - x^\ast) \\)（\\( J \\) の固有値 \\( \lambda \\) を指定）：単純混合は \\( w < 2/(1 - \lambda\_{\min}) \\) で収束、超えると発散 | 収束・発散の判定 |
| 線形模型 | 同じ模型で、履歴を十分長くした Pulay 法が、次元 + 数回の反復で収束すること | 反復回数 |
| 線形模型 | Pulay の係数：縁付き 1 次方程式の解と閉じた形 \\( A^{-1}\mathbf 1/(\mathbf 1^{\mathsf T}A^{-1}\mathbf 1) \\) の一致、\\( \sum\alpha = 1 \\) | \\( 10^{-8} \\) |
| 線形模型 | Thomas–Fermi 型応答 \\( \tilde R \propto -(q\_{\mathrm{TF}}^2/G^2) \\)：教科書の Kerker 混合（\\( q\_0 = q\_{\mathrm{TF}} \\)、\\( w = 1 \\)）が 1 反復で収束し、単純混合（\\( w = 0.3 \\)）が発散すること | 反復回数 |
| 解析 | 外挿係数：座標が時刻の 2 次式で動くとき、制約つき・制約なしの両方で \\( (3, -3, 1) \\) | \\( 10^{-8} \\) |
| 解析 | Kerker 因子の自動決定：立方格子で \\( q\_0 = 0.5 \\) bohr⁻¹、正方晶・斜方晶で上の式の値 | \\( 10^{-12} \\) |
| 照合 | 既存実装の fixture：小さな系（H₂O 分子、Si 結晶、Fe 結晶の共線、非共線）で、固有値・固有ベクトルを与えたときの \\( \mu \\)、\\( \rho(\mathbf L) \\)、\\( E(\mathbf L) \\) | \\( \mu \\)：\\( 10^{-8} \\) Hartree（ギャップ内の不定性を除く）、\\( \rho \\)・\\( E \\)：\\( 10^{-10} \\) |
| 照合 | 同じ fixture で、入力・出力の密度行列（または密度）の組と履歴を与えたときの、1 回のミキシングの結果と残差ノルム | \\( 10^{-10} \\)（相対） |
| 回帰 | 各方式で小さな系の SCF を収束させ、全エネルギーが方式によらず一致すること、反復回数が既存実装と同程度であること | 全エネルギー \\( 10^{-6} \\) Hartree、反復回数は目安 |

## 〔未確定〕

- `RMM-ADIIS`（拡張 Roothaan–Hall のエネルギー関数を使う DIIS）の具体的な式と、マニュアルでの扱い。
- O(N) 法（DC、Krylov）や NEGF での化学ポテンシャルの決め方と密度行列の作り方は [O(N) 法](on-methods.md)・[NEGF による電子輸送](negf.md) で扱いました（この章は Cluster と Band を対象としています）。
- 非共線の格子上の密度で、\\( \mathrm{Im}\\,n^{\uparrow\downarrow} \\) と磁化 \\( m\_y \\) の符号の対応（初期密度では \\( \mathrm{Im}\\,n^{\uparrow\downarrow} = -m\_y/2 \\)）。
  密度行列の定義からの導出と一致するかを、格子の章とあわせて確認する必要があります。
- 再開ファイルと外挿の履歴ファイルの形式（このノートで新しく定めてよい、〔設計〕に移す候補）。
- `empty.occupation.*`・`empty.states.*` など、特定の状態の占有を強制的に 0 にする入力（内殻正孔の計算などで使う）の扱い。
