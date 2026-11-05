# 二中心積分：重なり・運動エネルギー

> 状態：試作（レビュー待ち）。式は技術ノート 1（Eq. 36–46）の方法を、このノートの
> [Fourier 変換の規約](../design/notation.md) で導出し直したものです。

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

### 5. 特別な場合

- \\( \tau = 0 \\)（同じ原子上、\\( \mathbf L = 0 \\)）：\\( j\_L(0) = \delta\_{L0} \\) なので \\( L = 0 \\) の項だけが残り、
  \\( S\_{ab}(0) = \delta\_{ll'}\delta\_{mm'}\\,\frac{2}{\pi}\int\_0^\infty dk\\, k^2 \tilde R\_a \tilde R\_b = \delta\_{ll'}\delta\_{mm'}\int\_0^\infty dr\\, r^2 R\_a R\_b \\)。
  方向 \\( \hat{\boldsymbol\tau} \\) は不定なので、実装では特別に扱います。
- \\( \tau \ge r\_{c,a} + r\_{c,b} \\)：厳密に 0。波数空間の数値積分では小さな残差が出るので、
  この範囲の原子対は計算しません（疎ブロック形式の原子対の範囲と一致させます）。

### 6. 原子間ベクトルによる微分（力の計算用）

\\( f\_L(\tau) = I^{(n)}\_L(\tau) \\) と \\( Y^{\ast}\_{LM}(\hat{\boldsymbol\tau}) \\) の積の勾配は

$$
\nabla\_{\boldsymbol\tau}\left[f\_L(\tau)\\,Y^{\ast}\_{LM}(\hat{\boldsymbol\tau})\right]
= f\_L'(\tau)\\,Y^{\ast}\_{LM}\\,\hat{\boldsymbol\tau}
+ \frac{f\_L(\tau)}{\tau}\left(\partial\_\theta Y^{\ast}\_{LM}\\,\mathbf e\_\theta + \frac{1}{\sin\theta}\partial\_\varphi Y^{\ast}\_{LM}\\,\mathbf e\_\varphi\right)
$$

で、\\( f\_L'(\tau) = \int dk\\, k^{3+n} j\_L'(k\tau)\tilde R\_a\tilde R\_b \\) です。技術ノート 1 は、
力の計算でこれらの微分を解析的に評価すると述べています。極（\\( \sin\theta \to 0 \\)）では
\\( \frac{1}{\sin\theta}\partial\_\varphi Y \\) の極限値を使います。

## 数値パラメータ

〔公開〕技術ノート 15 は、波数空間の二中心積分に「カットオフエネルギー 3600 Ry、
Gauss–Legendre 積分点 128」を用いた例を示しています。3600 Ry は \\( k\_{\max} = 60 \\) bohr⁻¹ に当たります。
マニュアルには、既定値 3600 Ry の `1DFFT.EnergyCutoff` というキーワードがあります。

〔公開〕動径関数の球 Bessel 変換には、Toyoda–Ozaki の高速球 Bessel 変換（CPC **181**, 277）が使えます。

〔未確定〕

- 既存実装が実際に使う波数格子（点数、上限、積分則）と、`1DFFT.EnergyCutoff` との関係。
- \\( \tilde R(k) \\) を表として持ち補間するか、毎回計算するか（精度と速度の兼ね合いで〔設計〕として決めてよい）。

## テスト

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | 2 つの規格化 s 型 Gauss 関数 \\( e^{-\alpha r^2} \\)、\\( e^{-\beta r^2} \\)：\\( S = \left(\frac{2\sqrt{\alpha\beta}}{\alpha+\beta}\right)^{3/2} e^{-\mu\tau^2} \\)、\\( \mu = \frac{\alpha\beta}{\alpha+\beta} \\) | \\( 10^{-10} \\) |
| 解析 | 同じ組の運動エネルギー：\\( T = \mu\\,(3 - 2\mu\tau^2)\\, S \\) | \\( 10^{-10} \\) |
| 解析 | p・d 型 Gauss 関数の重なり（Gauss 積の定理による閉じた式） | \\( 10^{-10} \\) |
| 解析 | \\( \tau = 0 \\) での \\( S\_{ab} = \delta\_{ll'}\delta\_{mm'}\int r^2 R\_a R\_b\\,dr \\) | \\( 10^{-12} \\) |
| 不変量 | \\( S\_{ab}(\boldsymbol\tau) = \overline{S\_{ba}(-\boldsymbol\tau)} \\)、\\( T \\) も同様 | 丸め誤差 |
| 不変量 | 構造全体の回転で、行列が実球面調和関数の回転行列で変換されること | \\( 10^{-10} \\) |
| 不変量 | \\( \tau \ge r\_{c,a}+r\_{c,b} \\) で 0 | 厳密 |
| 不変量 | 解析微分と有限差分の一致 | \\( 10^{-6} \\) 程度 |
| 照合 | 実際の PAO（例：C6.0、O6.0、Fe6.0）で、いくつかの \\( \boldsymbol\tau \\) について fixture と比較 | 要決定 |
