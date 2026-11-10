# 二中心積分：重なり・運動エネルギー

> 状態：試作（レビュー待ち）。式は技術ノート 1（Eq. 36–46）の方法を、このノートの
> [Fourier 変換の規約](../../../conventions/notation.md) で導出し直したものです。

## 定義

原子 \\( I \\) の基底関数を \\( \phi\_a \\)、原子 \\( J \\) の基底関数を \\( \phi\_b \\) とし、
原子間ベクトルを \\( \boldsymbol\tau = \mathbf R\_J + \mathbf L - \mathbf R\_I \\) とします。

$$
S\_{ab}(\boldsymbol\tau) = \int d^3 r\\, \phi\_a^{\ast}(\mathbf r)\\, \phi\_b(\mathbf r - \boldsymbol\tau),\qquad
T\_{ab}(\boldsymbol\tau) = \int d^3 r\\, \phi\_a^{\ast}(\mathbf r) \left(-\tfrac12 \nabla^2\right) \phi\_b(\mathbf r - \boldsymbol\tau)
$$

## 方法〔公開〕

技術ノート 1 と同様に、波数空間で計算します。

### 1. 動径関数の変換

複素球面調和関数で \\( \phi\_a(\mathbf r) = R\_a(r) Y\_{lm}(\hat{\mathbf r}) \\) と書くと、平面波展開

$$
e^{-i\mathbf k\cdot\mathbf r} = 4\pi \sum\_{LM} (-i)^L j\_L(kr)\\, Y\_{LM}(\hat{\mathbf k})\\, Y^{\ast}\_{LM}(\hat{\mathbf r})
$$

と球面調和関数の正規直交性から

$$
\tilde\phi\_a(\mathbf k) = 4\pi (-i)^{l}\\, Y\_{lm}(\hat{\mathbf k})\\, \tilde R\_a(k),\qquad
\tilde R\_a(k) = \int\_0^\infty dr\\, r^2 j\_l(kr) R\_a(r)
$$

が得られます。

### 2. 重なり

Parseval の等式 \\( \int f^{\ast} g = (2\pi)^{-3}\int \tilde f^{\ast}\tilde g \\) と、平行移動
\\( \widetilde{\phi\_b(\cdot - \boldsymbol\tau)}(\mathbf k) = e^{-i\mathbf k\cdot\boldsymbol\tau}\tilde\phi\_b(\mathbf k) \\) を使い、
\\( e^{-i\mathbf k\cdot\boldsymbol\tau} \\) を再び平面波展開すると、

$$
S\_{ab}(\boldsymbol\tau) = 8 \sum\_{LM} i^{\\,l - l' - L}\\; G(lm;\\, l'm',\\, LM)\\; Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\\; I^{(0)}\_{L}(\tau)
$$

$$
I^{(n)}\_{L}(\tau) = \int\_0^\infty dk\\, k^{2+n}\\, j\_L(k\tau)\\, \tilde R\_a(k)\\, \tilde R\_b(k)
$$

となります。ここで \\( l' m' \\) は \\( \phi\_b \\) の量子数、\\( \tau = \lVert\boldsymbol\tau\rVert \\)、

$$
G(lm;\\, l'm',\\, LM) = \int d\Omega\\; Y^{\ast}\_{lm}\\, Y\_{l'm'}\\, Y\_{LM}
$$

は Gaunt 係数です。係数 8 は \\( (4\pi)^3/(2\pi)^3 \\) です。

- Gaunt 係数は \\( l + l' + L \\) が偶数かつ三角条件 \\( |l - l'| \le L \le l + l' \\) を満たすときだけ 0 でないので、
  \\( i^{\\,l-l'-L} \\) は \\( \pm 1 \\) です。和は有限です。
- \\( M = m - m' \\) 以外の項は 0 です。

### 3. 運動エネルギー

\\( -\tfrac12\nabla^2 \\) は波数空間で \\( k^2/2 \\) を掛けることに当たるので、

$$
T\_{ab}(\boldsymbol\tau) = 4 \sum\_{LM} i^{\\,l - l' - L}\\; G(lm;\\, l'm',\\, LM)\\; Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\\; I^{(2)}\_{L}(\tau)
$$

です。

> **検算**：重なりの式（\\( l, l' \le 2 \\) の複数の組）と運動エネルギーの式（s–s）は、
> Gauss 型動径関数 \\( r^l e^{-\alpha r^2} \\) について、三次元格子上の直接積分および解析解と
> 6 桁以上で一致することを `tools/check_two_center.py` で確認しています。

> **技術ノート 1 Eq. 46 についての注記**：同式では球 Bessel 関数の引数が \\( kr \\) と印字されていますが、
> 重なりの式（Eq. 44）との対応から \\( k\lvert\boldsymbol\tau\rvert \\) が正しいと考えられます。
> 上の導出もそれと整合します。

### 4. 実球面調和関数への変換

実基底 \\( \phi^{R}\_{\mu} = \sum\_m U^{(l)}\_{\mu m}\phi\_{m} \\)（[PAO 基底](pao-basis.md)）では

$$
S^{R}\_{\mu\nu} = \sum\_{m m'} \overline{U^{(l)}\_{\mu m}}\\; U^{(l')}\_{\nu m'}\\; S\_{m m'}
$$

です。結果は実数になります（虚部は丸め誤差のみ）。実 Gaunt 係数と実球面調和関数
\\( S\_{LM}(\hat{\boldsymbol\tau}) \\) を使って直接実数で計算しても構いません。

〔規約〕既存実装は複素基底で計算してから上の変換を行い、実部だけを残します。

## 原子対の範囲〔規約〕

重なりと運動エネルギーを計算する原子対は、

$$
\lVert \boldsymbol\tau \rVert \le r\_{c,I} + r\_{c,J}
$$

を満たす \\( (I, J, \mathbf L) \\) のすべてです。

- \\( r\_c \\) は各原子種の PAO ファイルの `radial.cutoff.pao` の値です（[PAO 基底](pao-basis.md)）。
- 不等号は等号を含みます。ちょうど \\( r\_{c,I} + r\_{c,J} \\) の距離の対も含めますが、値は 0 です。
- 原子自身（\\( J = I \\)、\\( \mathbf L = 0 \\)）も対に含めます。
- 周期像（\\( J = I \\)、\\( \mathbf L \ne 0 \\)）も、条件を満たせば別の対として扱います。
- 範囲外の対では \\( S = T = 0 \\) として行列に持ちません。

〔設計〕[疎ブロック形式](../design/sparse-block-format.md) の叩き台は狭義の不等号 \\( < \\) です。
境界上の対の値は 0 なので、どちらを採っても物理量は変わりません。
既存実装と原子対の並び・数まで照合する場合だけ、\\( \le \\) に合わせる必要があります。

## 波数積分の数値パラメータ〔規約〕

### 入力キーワード〔公開〕

マニュアル 8.2 節に次の 3 つのキーワードがあります。

| キーワード | 既定値 | 意味（マニュアルの説明を要約） |
|---|---|---|
| `1DFFT.EnergyCutoff` | 3600（Ry） | 動径関数の Fourier 変換で使う波数の上限をエネルギーで与える |
| `1DFFT.NumGridK` | 900 | 変換後の動径関数を保存する波数格子の点数 |
| `1DFFT.NumGridR` | 900 | 変換の数値積分に使う実空間の動径格子の点数 |

波数の上限は \\( k\_{\max} = \sqrt{E\_{\mathrm{cut}}[\mathrm{Ry}]} \\) です（[記法](../../../conventions/notation.md)）。既定値では 60 bohr⁻¹ です。

### 既存実装での使われ方

以下、\\( N\_K \\) を `1DFFT.NumGridK`、\\( N\_R \\) を `1DFFT.NumGridR` とします。

**(a) 動径関数の球 Bessel 変換** \\( \tilde R(k) = \int dr\\, r^2 j\_l(kr) R(r) \\)

| 項目 | 値 |
|---|---|
| 積分範囲 | \\( r\_{\min} \le r \le r\_c + 0.5 \\) bohr（\\( r\_{\min} \\) は PAO ファイルの最初の格子点） |
| 積分則 | \\( r \\) について等間隔の \\( N\_R \\) 区間、台形則。各点の \\( R(r) \\) は表の補間 |
| 前因子 | なし（[記法](../../../conventions/notation.md) の定義どおり。\\( \sqrt{2/\pi} \\) などは付けない） |
| 保存する波数 | \\( k\_j = j\\, k\_{\max}/N\_K \\)、\\( j = 0, 1, \dots, N\_K - 1 \\) |
| 保存の外 | \\( k > k\_{N\_K-1} \\) では \\( \tilde R = 0 \\)。間は補間 |
| 対象 | 基底指定で使う動径関数だけ（各 \\( l \\) の先頭 \\( N\_l \\) 個） |

**(b) 重なり・運動エネルギーの波数積分** \\( I^{(n)}\_L(\tau) = \int dk\\, k^{2+n} j\_L(k\tau)\tilde R\_a\tilde R\_b \\)

| 項目 | 値 |
|---|---|
| 積分範囲 | \\( k\_{\min} \le k \le k\_{\max} \\)、\\( k\_{\min} = 10^{-6} \\) bohr⁻¹ |
| 積分則 | \\( k \\) について等間隔の \\( N\_R \\) 区間（\\( N\_R + 1 \\) 点）、台形則 |
| \\( \tilde R(k) \\) | (a) の表の補間 |
| 球 Bessel 関数 | \\( L = 0, \dots, 2\max(l\_{\max,I}, l\_{\max,J}) \\) を各点で直接評価 |

- `1DFFT.NumGridR` は、(a) の \\( r \\) の区間数と (b) の \\( k \\) の区間数の両方に使われます。
  既定値では (b) の刻みは \\( 60/900 \approx 0.067 \\) bohr⁻¹ です。
- Gauss–Legendre 積分は使いません。技術ノート 15 の「Gauss–Legendre 128 点」は別の積分（射影子展開など）の話です。
- \\( k\_{\max} \\) 付近の \\( \tilde R \\) は十分小さいので、(a) の保存範囲が \\( k\_{\max} \\) に 1 刻み届かないことの影響は無視できます。

〔設計〕再実装では、同じ精度が出る限り積分則と格子は自由に選んでよいです
（例：Toyoda–Ozaki の高速球 Bessel 変換（CPC **181**, 277）、Gauss 型の積分則）。
ただし、上限 \\( k\_{\max} \\) は `1DFFT.EnergyCutoff` から上の式で決め、既定値 3600 Ry を守ります。
`1DFFT.NumGridK`、`1DFFT.NumGridR` は入力として受け付け、積分の分解能の指定と解釈します。
\\( \tilde R(k) \\) は原子種ごとに一度だけ表にして、原子対ごとには再計算しません。

**精度の目安**：実際の PAO（C6.0、H5.0、O6.0 の最初の s 関数）で、上の既定の格子による s–s 重なりと、
実空間での直接積分を比べると、\\( 0 \le \tau \le 12 \\) bohr で差は \\( 4\times 10^{-7} \\) 以下でした。
\\( \tau \ge r\_{c,a}+r\_{c,b} \\) での波数積分の残差は \\( 10^{-8} \\) 程度です（`tools/check_pao_file.py`）。

## 特別な場合

### \\( \tau = 0 \\)〔規約〕

\\( j\_L(0) = \delta\_{L0} \\) なので \\( L = 0 \\) の項だけが残り、

$$
S\_{ab}(0) = \delta\_{ll'}\delta\_{mm'}\\,\frac{2}{\pi}\int\_0^\infty dk\\, k^2 \tilde R\_a \tilde R\_b = \delta\_{ll'}\delta\_{mm'}\int\_0^\infty dr\\, r^2 R\_a R\_b,
\qquad
T\_{ab}(0) = \delta\_{ll'}\delta\_{mm'}\\,\frac{1}{\pi}\int\_0^\infty dk\\, k^4 \tilde R\_a \tilde R\_b
$$

です。方向 \\( \hat{\boldsymbol\tau} \\) は不定ですが、\\( L = 0 \\) の項は方向によりません。

- 既存実装は、\\( \tau < 10^{-10} \\) bohr を \\( \tau = 10^{-10} \\) bohr に置き換えて一般の式で計算します。
  このとき \\( L \ge 1 \\) の項は実質 0 になり、上の式と一致します。方向は任意の値で構いません。
- 原子自身の対（\\( J = I \\)、\\( \mathbf L = 0 \\)）では、原子間ベクトルによる微分（下記）を 0 とします。
- 〔設計〕再実装では \\( \tau < 10^{-10} \\) bohr を \\( \tau = 0 \\) とみなし、上の閉じた式を使います。
  同じ原子上で同じ \\( l \\) の異なる \\( \zeta \\) の重なりは、表の直交性（[PAO 基底](pao-basis.md)）から 0 に近い値になります。

### \\( \tau \ge r\_{c,a} + r\_{c,b} \\)

厳密に 0 です。波数空間の数値積分では小さな残差（上記）が出るので、この範囲の原子対は計算しません。

### 重なり行列に加える微小な乱数〔規約〕

既存実装は、重なり行列の各要素に、幅 \\( 10^{-13} \\) の一様乱数（\\( [-5\times10^{-14}, 5\times10^{-14}] \\)）を加えます。
数値的な縮退を崩すためのものと考えられます。

〔設計〕再実装では加えません。照合の許容誤差は \\( 10^{-12} \\) より大きくとります。

## 原子間ベクトルによる微分（力の計算用）

〔公開〕\\( f\_L(\tau) = I^{(n)}\_L(\tau) \\) と \\( Y^{\ast}\_{LM}(\hat{\boldsymbol\tau}) \\) の積の勾配は

$$
\nabla\_{\boldsymbol\tau}\left[f\_L(\tau)\\,Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\right]
{}= f\_L'(\tau)\\,Y^{\ast}\_{LM}\\,\hat{\boldsymbol\tau}
{}+ \frac{f\_L(\tau)}{\tau}\left(\partial\_\theta Y^{\ast}\_{LM}\\,\mathbf e\_\theta + \frac{1}{\sin\theta}\partial\_\varphi Y^{\ast}\_{LM}\\,\mathbf e\_\varphi\right)
$$

で、\\( f\_L'(\tau) = \int dk\\, k^{3+n} j\_L'(k\tau)\tilde R\_a\tilde R\_b \\) です。技術ノート 1 は、
力の計算でこれらの微分を解析的に評価すると述べています。

〔規約〕

- 既存実装は、\\( S \\)・\\( T \\) と同時に、中心原子の位置による微分
  \\( \partial S\_{ab}/\partial \mathbf R\_I = -\nabla\_{\boldsymbol\tau} S\_{ab} \\)（\\( T \\) も同様）の 3 成分を、同じ波数格子で計算して保持します。
- 極（\\( \lvert\sin\theta\rvert < 10^{-13} \\) 程度。\\( \boldsymbol\tau \\) が \\( z \\) 軸に平行な対）では、
  方位角を 0 とし、\\( \frac{1}{\sin\theta}\partial\_\varphi \\) の項を落とします。
  \\( Y\_{LM} \propto \sin^{\lvert M\rvert}\theta \\) なので、この項の極限値は \\( \lvert M\rvert = 1 \\) のときだけ 0 でなく、
  その寄与（方位角 0 では \\( y \\) 成分）が欠けます。例えば s 軌道と \\( p\_y \\) 軌道の対を \\( z \\) 軸上に置くと、
  正しい \\( \partial S/\partial \tau\_y = f(\tau)/\tau \\)（\\( f \\) は動径部分）は 0 でありませんが、この扱いでは 0 になります。
- 〔設計〕再実装では極限値を正しく入れるか、Cartesian 成分で微分を書いて特異点を避けます。
  そのため、既存実装との微分の照合では、\\( z \\) 軸に平行な原子対を含む構造を避けるか、この差を許します。
  〔未確定〕この差が既存実装の力にどれだけ現れるか（対称性で打ち消される場合が多い）は、
  [エネルギーと力](energy-force.md) の章で有限差分により確認します。
- 同じ処理の中で軌道角運動量演算子の行列要素も計算されますが、このノートの範囲外です。

## テスト

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 2 つの規格化 s 型 Gauss 関数 \\( e^{-\alpha r^2} \\)、\\( e^{-\beta r^2} \\)：\\( S = \left(\frac{2\sqrt{\alpha\beta}}{\alpha+\beta}\right)^{3/2} e^{-\mu\tau^2} \\)、\\( \mu = \frac{\alpha\beta}{\alpha+\beta} \\) | \\( 10^{-10} \\) |
| 解析 | 同じ組の運動エネルギー：\\( T = \mu\\,(3 - 2\mu\tau^2)\\, S \\) | \\( 10^{-10} \\) |
| 解析 | p・d 型 Gauss 関数の重なり（Gauss 積の定理による閉じた式） | \\( 10^{-10} \\) |
| 解析 | \\( \tau = 0 \\) での \\( S\_{ab} = \delta\_{ll'}\delta\_{mm'}\int r^2 R\_a R\_b\\,dr \\) | \\( 10^{-12} \\) |
| 解析 | \\( \tau < 10^{-10} \\) bohr の対で \\( \tau = 0 \\) の式と一致し、方向によらないこと | \\( 10^{-12} \\) |
| 不変量 | \\( S\_{ab}(\boldsymbol\tau) = \overline{S\_{ba}(-\boldsymbol\tau)} \\)、\\( T \\) も同様 | 丸め誤差 |
| 不変量 | 構造全体の回転で、行列が実球面調和関数の回転行列で変換されること | \\( 10^{-10} \\) |
| 不変量 | \\( \tau \ge r\_{c,a}+r\_{c,b} \\) の対が一覧にないこと。\\( \tau = r\_{c,a}+r\_{c,b} \\) ちょうどの対は既存実装では一覧にある | 厳密 |
| 不変量 | 解析微分と有限差分の一致（極の近くの \\( \boldsymbol\tau \\) を含む） | \\( 10^{-6} \\) 程度 |
| 収束 | 実際の PAO で、既定の波数格子と細かい格子（または実空間の直接積分）の差 | \\( 10^{-6} \\) |
| 照合 | 実際の PAO（例：C6.0、O6.0、Fe6.0S）で、いくつかの \\( \boldsymbol\tau \\) について既存実装の出力（scfout の重なり行列など）と比較 | \\( 10^{-6} \\)（積分則の差を許す。乱数の寄与 \\( 10^{-13} \\) は無視できる） |

> **検算**：解析解との照合は `tools/check_two_center.py`（上記）、実際の PAO での既定格子の精度は
> `tools/check_pao_file.py`（PAO ファイルのパスを引数に与える）で確認しています。
