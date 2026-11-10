# 非局所擬ポテンシャル

> 状態：試作（レビュー待ち）。分離型の式は原論文と技術ノート 2 を、このノートの記法で導出し直したものです。
> スピン軌道の行列要素と j 平均の規則は既存実装の挙動として確認し、`tools/check_nonlocal.py` で検算しています。

## 概要〔公開〕

- OpenMX は、Morrison–Bylander–Kleinman（MBK）形式の分離型ノルム保存擬ポテンシャルを使います
  （マニュアル 13.1 節。MBK は Vanderbilt の超ソフト擬ポテンシャルのノルム保存版と位置づけられています）。
  各 \\( l \\)（完全相対論では各 \\( j \\)）に複数の射影子を持てます。
- 擬ポテンシャルは、局所部分 \\( V\_{\mathrm{loc}}(r) \\) と分離型の非局所部分 \\( \hat V\_{\mathrm{NL}} \\) の和です。
  局所部分は [中性原子ポテンシャル](neutral-atom.md) に組み込まれるので、この章では非局所部分だけを扱います
  （技術ノート 1 Eq. 17–20, 23, 31）。
- 擬ポテンシャルのデータは原子種ごとに VPS ファイル（拡張子 `.vps`）で与えられます。
  原子種との対応は入力の `Definition.of.Atomic.Species` の第 3 列で指定します（マニュアル 13.1 節）。
  VPS ファイルは ADPACK が生成します（ADPACK マニュアル「Output files」）。

## 1. 分離型の形

### 1.1 Kleinman–Bylander と Blöchl〔公開〕

半局所型の擬ポテンシャル \\( V\_l(r) \\) から局所部分を引いた \\( \Delta V\_l = V\_l - V\_{\mathrm{loc}} \\) を、
参照状態 \\( \psi\_l \\) について再現する分離型にしたものが Kleinman–Bylander（KB）形式です
（Kleinman–Bylander, PRL **48**, 1425 (1982)）。

$$
\hat V\_{\mathrm{NL}} = \sum\_{lm} \frac{\lvert \Delta V\_l\\, \psi\_l\\, S\_{lm} \rangle \langle S\_{lm}\\, \psi\_l\\, \Delta V\_l \rvert}{\langle \psi\_l \rvert \Delta V\_l \lvert \psi\_l \rangle}
$$

Blöchl（PRB **41**, 5414 (1990)）は、1 つの \\( l \\) に複数の動径関数 \\( \bar R\_{l\zeta} \\) を使って精度を上げる形に拡張しました。
技術ノート 2 の Eq. 47–53 は、Gram–Schmidt で直交化した \\( \bar R\_{l\zeta} \\) と
\\( c\_{l\zeta} = \int r^2 dr\\, \bar R\_{l\zeta}\\, \Delta V\_l\\, \bar R\_{l\zeta} \\) を使い、
射影子 \\( \Delta V\_l \bar R\_{l\zeta} \\) と係数 \\( 1/c\_{l\zeta} \\) の対角和で書いています。
ADPACK では `Blochl.projector.num` が 1 つの \\( l \\) あたりの射影子数で、1 が KB 形式に当たります
（ADPACK マニュアル「Separable form」）。

### 1.2 MBK 形式〔公開〕

MBK 形式（Morrison–Bylander–Kleinman, PRB **47**, 6728 (1993)）では、同じ \\( l \\) の複数の参照エネルギー
\\( \varepsilon\_i \\) における擬波動関数 \\( \psi\_i \\) から

$$
\chi\_i = (\varepsilon\_i - \hat T - V\_{\mathrm{loc}})\\, \psi\_i,\qquad
B\_{ij} = \langle \psi\_i \vert \chi\_j \rangle
$$

を作り、

$$
\hat V\_{\mathrm{NL}} = \sum\_{ij} \lvert \chi\_i \rangle (B^{-1})\_{ij} \langle \chi\_j \rvert
$$

とします。一般化ノルム保存条件が成り立つと \\( B \\) は対称になります。参照エネルギーは ADPACK の入力
`pseudo.NandL` の第 5 列で与えます（ADPACK マニュアル「How the MBK scheme is different from the others」）。
1 つの \\( l \\) に状態が 1 つしかない場合は TM 形式に切り替わり、射影子数は `Blochl.projector.num` で決まります（同節）。

**対角化した形**〔設計〕：対称行列 \\( B \\) を \\( B = U \operatorname{diag}(b\_k) U^{\mathsf T} \\) と対角化し、

$$
\beta\_k = \frac{1}{\sqrt{\lvert b\_k \rvert}} \sum\_i \chi\_i U\_{ik},\qquad \varepsilon\_k = \operatorname{sign}(b\_k)
$$

と置くと \\( \hat V\_{\mathrm{NL}} = \sum\_k \lvert \beta\_k \rangle \varepsilon\_k \langle \beta\_k \rvert \\) になります
（`tools/check_nonlocal.py` で数値確認）。KB・Blöchl・MBK のどれも、この「射影子と対角係数」の形に書けます。

### 1.3 このノートでの表し方〔設計〕

原子 \\( K \\) の非局所部分を、射影子の番号 \\( p = 1, \dots, N\_p \\) について

$$
\hat V^{K}\_{\mathrm{NL}} = \sum\_{p} \varepsilon\_p \sum\_{\mu=1}^{2l\_p+1} \lvert \beta\_{p\mu}^{K} \rangle \langle \beta\_{p\mu}^{K} \rvert,\qquad
\beta\_{p\mu}^{K}(\mathbf r) = b\_p(\lvert \mathbf r - \mathbf R\_K \rvert)\\, S\_{l\_p \mu}\left(\widehat{\mathbf r - \mathbf R\_K}\right)
$$

と書きます。

- \\( l\_p \\) は射影子 \\( p \\) の角運動量、\\( b\_p(r) \\) は動径関数、\\( \varepsilon\_p \\) は射影子エネルギー（対角係数）です。
- \\( S\_{l\mu} \\) は [PAO 基底](pao-basis.md) の実球面調和関数で、成分の並びも同じです。
- 係数行列 \\( D \\) は射影子について対角（\\( D\_{pq} = \varepsilon\_p \delta\_{pq} \\)）です。
  非対角の \\( D \\) は VPS ファイルの形式に現れないので扱いません。
- 単位：\\( \varepsilon\_p\\, b\_p(r)\\, b\_p(r') \\) が Hartree·bohr⁻³ になるように選びます。
  \\( \varepsilon\_p = \pm 1 \\) のファイルでは大きさは \\( b\_p \\) が担い、\\( b\_p \\) の単位は \\( \mathrm{Hartree}^{1/2}\\,\mathrm{bohr}^{-3/2} \\) です。

〔規約〕DFT\_DATA19 の MBK 形式のファイルでは、射影子エネルギーはすべて \\( \pm 1 \\) です。つまり上の対角化と
規格化を済ませた射影子が格納されています。ただし、古い形式のファイル（例：`C02.vps`）や、
2 種類のファイルを混ぜた仮想結晶用のファイル（例：`C5Si5.vps`、エネルギーが \\( \pm 0.5 \\)）では
\\( \pm 1 \\) 以外の値も現れます。実装は \\( \varepsilon\_p \\) を一般の実数として扱ってください。

## 2. VPS ファイルの形式〔規約〕

### 2.1 全体の構造

- テキスト形式です。PAO ファイルと同じく、行単位の「キーワード 値」と、`<名前` で始まり `名前>` で終わるブロックからなります。
  `#` 以降はコメント、大文字・小文字は区別しません。論理値は `on`／`yes`／`true` などと `off`／`no`／`false` などです。
- ファイルの前半には、生成に使った ADPACK の入力ファイルがそのまま写されています（`vps.type`、`number.vps`、
  `<pseudo.NandL`、`Blochl.projector.num` など）。全電子計算の固有値・全エネルギーの記録も続きます。
  新形式では、これらは読み込みに使いません。
- 後半に、本体が読むデータのブロックがあります。新形式は `<project.energies`、`<Pseudo.Potentials`、
  `<density.PCC` の 3 つです。
- 同じキーワードが複数回現れた場合は最初の値を使います。

### 2.2 読み込む項目

| キーワード・ブロック | 意味 | 既定値 |
|---|---|---|
| `AtomSpecies` | 原子番号。仮想原子では小数（PAO ファイルと一致する必要がある） | 必須 |
| `total.electron` | 擬ポテンシャル生成時の全電子数 | `AtomSpecies` |
| `valence.electron` | 価電子数 | 0 |
| `grid.num.output` | 動径格子の点数 \\( N\_r \\) | 必須 |
| `j.dependent.pseudo.potentials` | j 依存（完全相対論）のデータか | off |
| `charge.pcc.calc` | 部分内殻補正の密度ブロックを含むか | off |
| `<project.energies` | 射影子の数・角運動量・エネルギー（新形式） | — |
| `<Pseudo.Potentials` | 動径格子、局所ポテンシャル、射影子の動径関数 | 必須 |
| `<density.PCC` | 部分内殻補正の密度（`charge.pcc.calc` が on のとき） | — |

イオンの有効電荷は \\( Z\_v = \\) `valence.electron` \\( + \\) `AtomSpecies` \\( - \\) `total.electron` です
（イオン状態で生成したファイルにも対応するため）。局所ポテンシャルの漸近形 \\( V\_{\mathrm{loc}}(r) \to -Z\_v/r \\)
と整合します（例：`C_PBE19.vps` の最後の格子点で \\( rV\_{\mathrm{loc}} = -4 \\)）。

PAO ファイルと VPS ファイルの `AtomSpecies` が異なると、エラーで停止します（空原子を除く）。

### 2.3 `<project.energies` ブロック（新形式）

1 行目は射影子の数 \\( N\_p \\)、続く \\( N\_p \\) 行が各射影子 \\( p \\) の記述です。

- j 依存でないとき：各行は「\\( l\_p \\)　\\( \varepsilon\_p \\)」の 2 列。
- j 依存のとき：各行は「\\( l\_p \\)　\\( \varepsilon\_{p+} \\)　\\( \varepsilon\_{p-} \\)」の 3 列。
  \\( + \\) は \\( j = l\_p + 1/2 \\)、\\( - \\) は \\( j = l\_p - 1/2 \\) です。

エネルギーの単位は Hartree です。行の順序が射影子の番号 \\( p \\) を決めます。DFT\_DATA19 のファイルでは
\\( l \\) の昇順に、同じ \\( l \\) の射影子がまとまって並んでいますが、読み込み側は各行の \\( l\_p \\) を使い、
並び方を仮定しないでください。\\( N\_p = 0 \\) も許されます（空原子用の `E.vps`）。

### 2.4 `<Pseudo.Potentials` ブロック

\\( N\_r \\) 行からなり、各行は

| 列 | 内容 |
|---|---|
| 1 | \\( x = \ln r \\) |
| 2 | \\( r \\)（bohr） |
| 3 | 局所ポテンシャル \\( V\_{\mathrm{loc}}(r) \\)（Hartree） |
| 4 以降（j 依存でない） | \\( b\_1(r), b\_2(r), \dots, b\_{N\_p}(r) \\) |
| 4 以降（j 依存） | \\( b\_{1+}(r), b\_{1-}(r), b\_{2+}(r), b\_{2-}(r), \dots \\)（射影子ごとに \\( j = l+1/2 \\)、\\( j = l-1/2 \\) の順） |

です。列数は j 依存でなければ \\( 3 + N\_p \\)、j 依存なら \\( 3 + 2N\_p \\) です
（例：`C_PBE19.vps` は \\( N\_p = 12 \\)、j 依存で 27 列）。

- 格子は \\( \ln r \\) について等間隔です（`C_PBE19.vps` では \\( x \\) が −7.5 から 2.8 まで 500 点）。
- 射影子の値は \\( r\\, b\_p(r) \\) ではなく \\( b\_p(r) \\) そのものです。原点付近で \\( b\_p \propto r^{l\_p} \\) のように振る舞います。
  重なりは \\( \int\_0^\infty r^2 dr\\, b\_p(r) R(r) \\) の形で使います。
- 射影子は生成時のカットオフ半径の外で急速に 0 に近づきますが、表の最後まで値が格納されています
  （`C_PBE19.vps` では、最大値の \\( 10^{-6} \\) 倍を下回るのはおよそ 1.3〜3 bohr より外）。
  射影子の打ち切り半径はファイルに明示されていません。

### 2.5 `<density.PCC` ブロック

`charge.pcc.calc` が on のときだけ存在します。各行は「\\( \ln r \\)　\\( r \\)　\\( \rho\_{\mathrm{pcc}}(r) \\)」で、
格子は `<Pseudo.Potentials` と同じでなければなりません（一致しなければエラー）。
\\( \rho\_{\mathrm{pcc}} \\) は電子数密度（bohr⁻³）で、\\( \int 4\pi r^2 \rho\_{\mathrm{pcc}}\\, dr \\) が部分内殻電荷になります
（`C_PBE19.vps` で約 0.25）。この密度は交換相関の計算で価電子密度に加えるもので、
[実空間グリッド・Poisson・交換相関](grid-poisson-xc.md) で使います。off のときは 0 として扱います。

### 2.6 旧形式

`<project.energies` ブロックがないファイルは旧形式として読みます。DFT\_DATA19 には旧形式のファイルはありません。

- 旧形式は j 依存のデータを持てません。`j.dependent.pseudo.potentials` が on ならエラーです。
- `number.vps` が記述行の数、`<pseudo.NandL` の各行が「番号　主量子数　\\( l \\)　カットオフ半径」です。
- `local.part.vps` は、局所部分に使った行の番号（0 始まり）です。この行を除いた残りの行が、その順で
  1 つずつ射影子になります。番号が範囲外のときは、どの行も除きません。
- `<projection.energies` ブロックの各行が「番号　\\( \varepsilon\_p \\)」で、射影子と同じ順に並びます。
- `<Pseudo.Potentials` の列は新形式の j 依存でない場合と同じです。

### 2.7 補間と範囲外の値

- 補間の方式は原則として自由です（[PAO 基底](pao-basis.md) と同じく、ビット一致は求めません）。
- 格子の最大半径より外では 0 とします。
- 格子の最小半径より内側では、\\( l = 0 \\) の射影子は最初の区間からの外挿、\\( l \ge 1 \\) では 0 とします。
  最小半径は \\( e^{-7.5} \approx 5.5\times 10^{-4} \\) bohr 程度なので、結果への影響は無視できます。

## 3. 行列要素の組み立て

### 3.1 射影子との重なり〔公開〕

原子 \\( I \\) 上の基底関数 \\( \phi\_a \\)（\\( l\_a \\)、動径関数 \\( R\_a \\)）と、原子 \\( K \\) の射影子 \\( \beta^{K}\_{p\mu} \\) の重なりを

$$
C^{K}\_{a,p\mu}(\boldsymbol\tau) = \int d^3 r\\, \phi\_a(\mathbf r)\\, \beta\_{p\mu}(\mathbf r - \boldsymbol\tau),\qquad
\boldsymbol\tau = \mathbf R\_K + \mathbf L' - \mathbf R\_I
$$

と書きます。これは [二中心積分](two-center-overlap-kinetic.md) の重なり \\( S\_{ab}(\boldsymbol\tau) \\) で、
相手の動径関数 \\( R\_b \\) を \\( b\_p \\) に置き換えたものと同じです（同章の 2 節の式と実球面調和関数への変換）。
基底関数も射影子も実関数なので \\( C \\) は実数です。

### 3.2 行列要素（スカラーの場合）〔公開〕

原子 \\( I \\) 上の \\( \phi\_a \\) と、格子ベクトル \\( \mathbf L \\) だけずれた原子 \\( J \\) 上の \\( \phi\_b \\) の間の行列要素は

$$
V^{\mathrm{NL}}\_{ab}(\mathbf L) = \sum\_{K,\mathbf L'} \sum\_{p} \varepsilon\_p \sum\_{\mu}
C^{K}\_{a,p\mu}(\mathbf R\_K + \mathbf L' - \mathbf R\_I)\\;
C^{K}\_{b,p\mu}(\mathbf R\_K + \mathbf L' - \mathbf R\_J - \mathbf L)
$$

です（技術ノート 1 Eq. 11, 78 の非局所項を、分離型を代入して書き下したもの）。
和は単位胞内のすべての原子 \\( K \\) とその周期像 \\( \mathbf L' \\) についてとります。
\\( K = I \\) や \\( K = J \\)（同じ周期像）の項も含みます。
スピンを考えない計算やスピン共線の計算では、この行列が上向き・下向きの両方に同じだけ加わります。

### 3.3 第三の原子の範囲〔規約〕

既存実装は、行列要素に寄与する第三の原子を、PAO のカットオフ半径（基底指定の `C6.0-…` の 6.0、
[PAO 基底](pao-basis.md)）だけで選びます。

- 原子対 \\( (I, J+\mathbf L) \\) は、重なり行列と同じ範囲 \\( \lvert \boldsymbol\tau\_{IJ} \rvert \le r\_{c,I} + r\_{c,J} \\) のものだけを計算します。
  この範囲の外の原子対の非局所行列要素は捨てます。
- 第三の原子 \\( K + \mathbf L' \\) は、\\( \lvert \boldsymbol\tau\_{IK} \rvert \le r\_{c,I} + r\_{c,K} \\) と
  \\( \lvert \boldsymbol\tau\_{JK} \rvert \le r\_{c,J} + r\_{c,K} \\) の両方を満たすものだけを和に入れます。
  ここで \\( r\_{c,K} \\) は原子 \\( K \\) の **PAO** のカットオフ半径で、射影子の広がりではありません。

射影子の広がり \\( r\_\beta \\) は PAO のカットオフ半径よりずっと小さいので、2 つ目の条件で落ちる項は実質的に 0 です。
一方、1 つ目の条件は近似です。物理的には \\( \lvert \boldsymbol\tau\_{IJ} \rvert \\) が \\( r\_{c,I} + r\_{c,J} + 2r\_\beta \\) まで
0 でない要素がありえますが、疎行列の形を重なり行列と揃えるために捨てています。
その切断を滑らかにするのが次の減衰関数です。

### 3.4 減衰関数〔規約〕

原子対 \\( (I, J+\mathbf L) \\) の非局所行列要素全体に、原子間距離 \\( d = \lvert \boldsymbol\tau\_{IJ} \rvert \\) の関数

$$
f(d) = \begin{cases}
1 & d < r\_{\mathrm{cut}} - 1 \\\\
10t^3 - 15t^4 + 6t^5,\quad t = r\_{\mathrm{cut}} - d & r\_{\mathrm{cut}} - 1 \le d \le r\_{\mathrm{cut}} \\\\
0 & d > r\_{\mathrm{cut}}
\end{cases},\qquad r\_{\mathrm{cut}} = r\_{c,I} + r\_{c,J}
$$

を掛けます（長さは bohr、緩和幅は 1 bohr）。\\( f \\) は両端で 2 階微分まで連続です（`tools/check_nonlocal.py`）。
同じ減衰関数は、中性原子ポテンシャルの射影子展開による行列要素にも使われています。

〔設計〕再実装はこの切断と減衰関数をそのまま採用します。既存実装の fixture と比較できること、
疎ブロック形式（[H(R)・S(R) 疎ブロック形式](../design/sparse-block-format.md)）の原子対の集合を
重なり行列と共有できることが理由です。切断しない厳密な計算は、検証用の別経路として持つことにします。

### 3.5 周期像

\\( \mathbf L \\)、\\( \mathbf L' \\) は、上の距離の条件を満たすすべての格子ベクトルをとります。
\\( V^{\mathrm{NL}}\_{ab}(\mathbf L) \\) は \\( H(\mathbf R) \\) の 1 ブロックとして蓄えます。
Hermite 性 \\( V^{\mathrm{NL}}\_{ab}(\mathbf L) = V^{\mathrm{NL}}\_{ba}(-\mathbf L) \\)（スカラーの場合は実対称）が成り立ちます。

## 4. スカラー相対論での j 依存データの扱い

### 4.1 方法〔公開〕

`scf.SpinOrbit.Coupling` が off（既定値）のとき、j 依存の擬ポテンシャルは j の縮退度を重みにして平均され、
スカラー相対論の擬ポテンシャルとして使われます（マニュアル 38.3 節）。
j の縮退度は \\( 2j+1 \\) なので、角運動量 \\( l \\) で \\( j = l \pm 1/2 \\) の重みは

$$
w\_+ = \frac{l+1}{2l+1},\qquad w\_- = \frac{l}{2l+1}
$$

です。マニュアルには、何を平均するか（演算子か、射影子とエネルギーか）は書かれていません。

### 4.2 既存実装の平均の仕方〔規約〕

既存実装は、射影子の動径関数とエネルギーを **別々に** 平均します。

$$
\bar b\_p(r) = w\_+ b\_{p+}(r) + w\_- b\_{p-}(r),\qquad
\bar\varepsilon\_p = w\_+ \varepsilon\_{p+} + w\_- \varepsilon\_{p-}
$$

そして \\( \bar\varepsilon\_p\\, \lvert \bar b\_p S \rangle\langle \bar b\_p S \rvert \\) を 3.2 節の式で使います。
共線スピン計算でも、スピン軌道相互作用なしの非共線計算でも同じです。

〔公開〕（導出）一方、j 依存の演算子（5 節）をスピンについて平均したもの（スピンのトレースの 1/2）は

$$
\sum\_p \left( w\_+ \varepsilon\_{p+} \lvert b\_{p+} S \rangle\langle b\_{p+} S \rvert + w\_- \varepsilon\_{p-} \lvert b\_{p-} S \rangle\langle b\_{p-} S \rvert \right)
$$

です（`tools/check_nonlocal.py` で \\( \tfrac12 \mathrm{Tr}\_\sigma P\_\pm = w\_\pm \\) を確認）。両者は \\( b\_{p+} = b\_{p-} \\) かつ
\\( \varepsilon\_{p+} = \varepsilon\_{p-} \\) なら一致しますが、一般には異なります。
DFT\_DATA19 には、同じ射影子で \\( \varepsilon\_{p+} = 1 \\)、\\( \varepsilon\_{p-} = -1 \\) のように符号が異なるもの
（例：`Ba_PBE19.vps`、`Ag_PBE19.vps` の \\( l = 1 \\)）があり、このとき既存実装の \\( \bar\varepsilon\_p \\) は \\( 1/(2l+1) \\) になります。

〔設計〕既定では既存実装と同じ平均（射影子とエネルギーを別々に平均）を使い、fixture と比較します。
演算子の平均は選択肢として用意し、両者の差（全エネルギー、バンド）を記録しておきます。

## 5. 完全相対論とスピン軌道相互作用

### 5.1 j 依存の非局所演算子〔公開〕

j 依存の擬ポテンシャルは、角運動量 \\( l \\) とスピン 1/2 を合成した全角運動量 \\( j = l \pm 1/2 \\) の
部分空間への射影を使って書かれます（Bachelet–Schlüter, PRB **25**, 2103 (1982)；
Theurich–Hill, PRB **64**, 073106 (2001)；技術ノート 2 Eq. 40–58）。
技術ノート 2 は Clebsch–Gordan 係数で作ったスピノル \\( \Phi\_{JM} \\) への射影の和で書いています。
\\( \mathbf L\cdot\boldsymbol\sigma = 2\\,\mathbf L\cdot\mathbf S \\) の固有値が \\( j = l+1/2 \\) で \\( l \\)、\\( j = l-1/2 \\) で \\( -(l+1) \\) であることから、
この射影は

$$
P^{(l)}\_{+} = \frac{(l+1)\\,\mathbb 1 + \mathbf L\cdot\boldsymbol\sigma}{2l+1},\qquad
P^{(l)}\_{-} = \frac{l\\,\mathbb 1 - \mathbf L\cdot\boldsymbol\sigma}{2l+1}
$$

と書けます（\\( P\_+ + P\_- = \mathbb 1 \\)、\\( P\_\pm^2 = P\_\pm \\)、\\( \mathrm{Tr}\\,P\_\pm = 2j+1 \\)）。
これを使うと原子 \\( K \\) の非局所演算子は

$$
\hat V^{K}\_{\mathrm{NL}} = \sum\_p \sum\_{s=\pm} \varepsilon\_{ps} \sum\_{\mu\nu} \sum\_{\sigma\sigma'}
\lvert \beta^{K}\_{ps\mu}\\, \sigma \rangle \left(P^{(l\_p)}\_{s}\right)^{\sigma\sigma'}\_{\mu\nu} \langle \beta^{K}\_{ps\nu}\\, \sigma' \rvert
$$

です。\\( \beta^{K}\_{ps\mu} \\) の動径関数は \\( b\_{ps} \\) です。\\( P \\) の行列要素は、実球面調和関数とスピンの直積基底
\\( \lvert S\_{l\mu} \rangle \lvert \sigma \rangle \\) でとります（\\( \sigma = \uparrow, \downarrow \\) は \\( \sigma\_z \\) の固有状態）。
上の \\( P\_\pm \\) が Clebsch–Gordan スピノルによる射影と一致することは、\\( l = 1, 2, 3 \\) について
`tools/check_nonlocal.py` で確認しています。

### 5.2 実球面調和関数での角運動量行列〔設計〕

\\( (L\_k)\_{\mu\nu} = \langle S\_{l\mu} \vert \hat L\_k \vert S\_{l\nu} \rangle \\)、\\( \hat{\mathbf L} = -i\\,\mathbf r\times\nabla \\) とし、
\\( L\_\pm = L\_x \pm i L\_y \\) とします。スピン空間で

$$
\mathbf L\cdot\boldsymbol\sigma = \begin{pmatrix} L\_z & L\_- \\\\ L\_+ & -L\_z \end{pmatrix}
$$

です（行・列ともに \\( \uparrow, \downarrow \\) の順）。実基底では \\( L\_z \\) は純虚数の反対称行列、\\( L\_- \\) は反対称行列で、
\\( L\_+ = L\_-^{\dagger} \\) です。[PAO 基底](pao-basis.md) の並びと符号で 0 でない上三角の要素は次のとおりです
（成分番号は 0 始まり。`tools/check_nonlocal.py` で計算）。

- \\( l = 1 \\)（\\( x, y, z \\)）：\\( (L\_z)\_{01} = -i \\)；\\( (L\_-)\_{02} = 1 \\)、\\( (L\_-)\_{12} = -i \\)。
- \\( l = 2 \\)（\\( 3z^2-r^2, x^2-y^2, xy, xz, yz \\)）：\\( (L\_z)\_{12} = -2i \\)、\\( (L\_z)\_{34} = -i \\)；
  \\( (L\_-)\_{03} = -\sqrt3 \\)、\\( (L\_-)\_{04} = \sqrt3\\, i \\)、\\( (L\_-)\_{13} = 1 \\)、\\( (L\_-)\_{14} = i \\)、
  \\( (L\_-)\_{23} = -i \\)、\\( (L\_-)\_{24} = 1 \\)。
- \\( l = 3 \\)（PAO 基底の表の順）：\\( (L\_z)\_{12} = -i \\)、\\( (L\_z)\_{34} = -2i \\)、\\( (L\_z)\_{56} = -3i \\)；
  \\( (L\_-)\_{01} = -\sqrt6 \\)、\\( (L\_-)\_{02} = \sqrt6\\, i \\)、
  \\( (L\_-)\_{13} = -\sqrt{5/2} \\)、\\( (L\_-)\_{14} = \sqrt{5/2}\\, i \\)、
  \\( (L\_-)\_{23} = -\sqrt{5/2}\\, i \\)、\\( (L\_-)\_{24} = -\sqrt{5/2} \\)、
  \\( (L\_-)\_{35} = -\sqrt{3/2} \\)、\\( (L\_-)\_{36} = \sqrt{3/2}\\, i \\)、
  \\( (L\_-)\_{45} = -\sqrt{3/2}\\, i \\)、\\( (L\_-)\_{46} = -\sqrt{3/2} \\)。

下三角は反対称性から決まります（\\( (L\_z)\_{\nu\mu} = -(L\_z)\_{\mu\nu} \\)、\\( (L\_-)\_{\nu\mu} = -(L\_-)\_{\mu\nu} \\)）。

〔規約〕既存実装のスピン軌道の行列要素は、この \\( L\_z \\)、\\( L\_- \\) と 5.1 節の \\( P\_\pm \\) で書いたものと
\\( l = 1, 2, 3 \\) のすべての要素で一致します（符号を含む）。

### 5.3 行列要素〔公開〕（導出）

3.2 節と同じ記号で、スピンの組 \\( (\sigma, \sigma') \\) の行列要素は

$$
V^{\mathrm{NL},\sigma\sigma'}\_{ab}(\mathbf L) = \sum\_{K,\mathbf L'} \sum\_{p} \sum\_{s=\pm} \varepsilon\_{ps} \sum\_{\mu\nu}
C^{K}\_{a,ps\mu}\\, \left(P^{(l\_p)}\_{s}\right)^{\sigma\sigma'}\_{\mu\nu}\\, C^{K}\_{b,ps\nu}
$$

です。ブロックごとに書くと（\\( C\_{a,s} \\) を \\( \mu \\) についての行ベクトルとして）

$$
V^{\uparrow\uparrow} = \sum \frac{\varepsilon\_{+}}{2l+1} C\_{a,+}\left[(l+1) + L\_z\right] C\_{b,+}^{\mathsf T}
{}+ \frac{\varepsilon\_{-}}{2l+1} C\_{a,-}\left[l - L\_z\right] C\_{b,-}^{\mathsf T}
$$

$$
V^{\downarrow\downarrow} = \sum \frac{\varepsilon\_{+}}{2l+1} C\_{a,+}\left[(l+1) - L\_z\right] C\_{b,+}^{\mathsf T}
{}+ \frac{\varepsilon\_{-}}{2l+1} C\_{a,-}\left[l + L\_z\right] C\_{b,-}^{\mathsf T}
$$

$$
V^{\uparrow\downarrow} = \sum \frac{\varepsilon\_{+}}{2l+1} C\_{a,+} L\_- C\_{b,+}^{\mathsf T}
{}- \frac{\varepsilon\_{-}}{2l+1} C\_{a,-} L\_- C\_{b,-}^{\mathsf T},\qquad
V^{\downarrow\uparrow}\_{ab}(\mathbf L) = \overline{V^{\uparrow\downarrow}\_{ba}(-\mathbf L)}
$$

です（和は \\( K, \mathbf L', p \\) について。3.3・3.4 節の範囲と減衰関数もそのまま適用）。

- 対角ブロックの実部は j 依存の重み付きの和、虚部が \\( L\_z \\) の項（スピン軌道の対角部分）です。
- \\( V^{\uparrow\downarrow} \\) は複素数です。既存実装は \\( \uparrow\uparrow \\)、\\( \downarrow\downarrow \\)、\\( \uparrow\downarrow \\) の 3 ブロックの
  実部と虚部を持ち、\\( \downarrow\uparrow \\) は Hermite 性から作ります〔規約〕。
- j 依存でない原子種（スカラーのファイル）が混ざる場合、その原子 \\( K \\) の寄与は 3.2 節の式で、
  \\( \uparrow\uparrow \\) と \\( \downarrow\downarrow \\) に同じだけ加わります〔規約〕。
- 既存実装は j 依存の射影子を \\( l \le 3 \\) についてだけ扱い、\\( l \ge 4 \\) があるとエラーで停止します〔規約〕。

### 5.4 入力キーワードと条件〔公開〕〔規約〕

- `scf.SpinOrbit.Coupling`（on/off、既定 off）：スピン軌道相互作用を自己無撞着に入れます（マニュアル 38.1 節）。
  j 依存の VPS ファイルを使います。
- 〔規約〕on にできるのは非共線計算（`scf.SpinPolarization` が `NC`）だけです。
  共線のスピン分極計算やスピンなしの計算で on にすると、入力エラーで停止します。
- 〔規約〕on でも、使う VPS ファイルがすべて j 依存でなければ、off に切り替えて計算を続けます。
- 〔規約〕off のとき、j 依存のファイルは 4 節の規則で平均されます。

### 5.5 スピン軌道相互作用の強さの調整〔公開〕〔規約〕

〔公開〕`<scf.SO.factor` ブロックで、原子種と \\( l \\)（s, p, d, f）ごとにスピン軌道相互作用の強さを変えられます
（マニュアル 38.2 節）。各行は「原子種名　s 係数　p 係数　d 係数　f 係数」の形で、各係数の前に
`s`、`p`、`d`、`f` の記号が置かれます。係数 1.0 が実際の原子の強さです。すべての原子種について f まで書く必要があります。
なお、ADPACK にも擬ポテンシャル生成時に強さを変える `<SO.factor` ブロックがありますが、これは別物です
（ADPACK マニュアル「Enhancement or depletion of a spin-orbit coupling」）。

〔規約〕

- このブロックは非共線計算で `scf.SpinOrbit.Coupling` が on のときだけ読まれます。ブロックがなければ係数はすべて 1 です。
- 各行の記号（`s` など）は読み飛ばされ、係数は \\( l = 0, 1, 2, 3 \\) の順に対応します。
- 係数 \\( \lambda\_l \\) は射影子の動径関数に作用します。\\( \bar b = (b\_+ + b\_-)/2 \\)、\\( \delta b = (b\_+ - b\_-)/2 \\) として

$$
b\_{p+} \to \bar b\_p + \lambda\_{l\_p}\\, \delta b\_p,\qquad b\_{p-} \to \bar b\_p - \lambda\_{l\_p}\\, \delta b\_p
$$

  と置き換え、エネルギー \\( \varepsilon\_{p\pm} \\) は変えません。

〔公開〕（導出）この方法では \\( \lambda = 0 \\) でもスピン軌道相互作用は完全には消えません。\\( b\_{p+} = b\_{p-} = \bar b\_p \\) のとき

$$
\varepsilon\_{+} P\_+ + \varepsilon\_{-} P\_- = \frac{(l+1)\varepsilon\_+ + l\varepsilon\_-}{2l+1}\\,\mathbb 1 + \frac{\varepsilon\_+ - \varepsilon\_-}{2l+1}\\,\mathbf L\cdot\boldsymbol\sigma
$$

なので、\\( \varepsilon\_{p+} \ne \varepsilon\_{p-} \\) の射影子では \\( \mathbf L\cdot\boldsymbol\sigma \\) の項が残ります（`tools/check_nonlocal.py`）。
再実装では既存実装と同じ置き換えを使い、この性質をテストで確認します。

## 6. 原子位置による微分（力）

〔公開〕非局所項のエネルギーは \\( E\_{\mathrm{NL}} = \sum \rho\\, V^{\mathrm{NL}} \\)（密度行列との縮約）で、
力への寄与は行列要素の微分と密度行列の積、および重なり行列の微分と エネルギー密度行列の積（Pulay 項。
[全エネルギーと力](energy-force.md) で扱う）に分かれます（技術ノート 1 Eq. 62, 71–75）。
技術ノート 1 は行列要素の微分を解析的に評価するとしています。ここでは前者を式にします。

〔公開〕（導出）非共線の場合も含めて

$$
E\_{\mathrm{NL}} = \sum\_{\mathbf L} \sum\_{ab} \sum\_{\sigma\sigma'} \rho^{\sigma\sigma'}\_{ab}(\mathbf L)\\, V^{\mathrm{NL},\sigma\sigma'}\_{ab}(\mathbf L)
$$

と書きます。複素共役は付けない積の和で、実部と虚部に分けると \\( \mathrm{Re}\\,\rho\\,\mathrm{Re}\\,V - \mathrm{Im}\\,\rho\\,\mathrm{Im}\\,V \\)、
スピン非対角 \\( \uparrow\downarrow \\) には係数 2 が付きます
（[非共線・スピン軌道計算のエネルギー・力・応力](noncollinear-energy-force.md) 1.4 節の縮約の規則。バンドエネルギーとの一致を数値で確認済みで、既存実装もこの縮約を使います〔規約〕）。
\\( \rho \\) の定義と位相の約束は [k 点と Bloch 和](kpoints-bloch.md)・[密度行列と電荷ミキシング](density-mixing.md) と [H(R)・S(R) 疎ブロック形式](../design/sparse-block-format.md) に従います。
（以前の版の \\( \rho^{\sigma'\sigma}\_{ba}(-\mathbf L) \\) との縮約は、これらの章の \\( \rho \\) の規約と虚部どうしの積の符号が合わないため、上の形に改めました。）
3.4 節の減衰関数を含めて \\( V\_{ab} = f(d\_{IJ})\\, h\_{ab} \\)（\\( h \\) は減衰前の和）とすると、原子 \\( X \\) の位置による微分は

$$
\frac{\partial V\_{ab}}{\partial \mathbf R\_X} = f'(d\_{IJ})\\, \frac{\partial d\_{IJ}}{\partial \mathbf R\_X}\\, h\_{ab} + f(d\_{IJ})\\, \frac{\partial h\_{ab}}{\partial \mathbf R\_X},\qquad
\frac{\partial d\_{IJ}}{\partial \mathbf R\_J} = -\frac{\partial d\_{IJ}}{\partial \mathbf R\_I} = \hat{\boldsymbol\tau}\_{IJ}
$$

です（\\( f' \\) は \\( d \\) についての微分。\\( t = r\_{\mathrm{cut}} - d \\) なので \\( f'(d) = -(30t^2 - 60t^3 + 30t^4) \\)）。
\\( h \\) の微分は、重なり \\( C^{K}\_{a}(\boldsymbol\tau) \\) が \\( \boldsymbol\tau = \mathbf R\_K + \mathbf L' - \mathbf R\_I \\) だけに依存することから

$$
\frac{\partial h\_{ab}}{\partial \mathbf R\_X} = \sum\_{K,\mathbf L'}\sum\_{p,s} \varepsilon\_{ps}
\left[ (\delta\_{XK} - \delta\_{XI})\\, \nabla\_{\boldsymbol\tau} C^{K}\_{a,ps}\\; P\\; C^{K\\,\mathsf T}\_{b,ps}
{}+ (\delta\_{XK} - \delta\_{XJ})\\, C^{K}\_{a,ps}\\; P\\; \nabla\_{\boldsymbol\tau} C^{K\\,\mathsf T}\_{b,ps} \right]
$$

です（スカラーの場合は \\( P \\) を単位行列、\\( s \\) の和を 1 項にする）。
\\( \nabla\_{\boldsymbol\tau} C \\) は [二中心積分](two-center-overlap-kinetic.md) の 6 節の式で計算します。

- 第三の原子 \\( K \\) 自身に力がかかる（\\( X = K \\) の項）ことに注意してください。
- \\( K = I \\) の項では \\( \delta\_{XK} - \delta\_{XI} = 0 \\) となり、同じ原子上の重なりは位置によらないことと整合します。
- 力の寄与は \\( \mathbf F\_X = -\partial E\_{\mathrm{NL}}/\partial \mathbf R\_X \\) です。系全体の並進で \\( \sum\_X \mathbf F\_X = 0 \\) になります。
- 〔規約〕既存実装の力の計算も減衰関数の微分の項を含みます。

## 7. 数値パラメータ〔規約〕

- \\( C^{K}\_{a,p\mu} \\) は二中心積分と同じく波数空間で計算します。動径関数の球 Bessel 変換は、
  実空間で最小半径から「原子 \\( K \\) の PAO カットオフ半径 + 0.5 bohr」までを積分範囲とし、
  `1DFFT.NumGridR`（既定 900）点の等間隔格子を使います。変換した関数は、0 から \\( k\_{\max} \\) までの
  `1DFFT.NumGridK`（既定 900）点の等間隔格子の表として持ちます。
- \\( k\_{\max} = \sqrt{E\_{\mathrm{cut}}} \\) bohr⁻¹ で、\\( E\_{\mathrm{cut}} \\) は `1DFFT.EnergyCutoff`（Ry、既定 3600）です。
- 波数についての積分（二中心積分の \\( I^{(0)}\_L \\) に当たるもの）は、\\( 10^{-6} \\) から \\( k\_{\max} \\) までを
  `1DFFT.NumGridR` 等分した格子上の台形則で行います。

ビット一致は求めないので、これらは〔設計〕で変えてよい既定値として扱い、fixture との差が許容誤差内に収まることを確認します。

## 〔未確定〕

- 射影子の打ち切り半径 \\( r\_\beta \\) を明示的に定めるか（第三の原子の選び方を PAO カットオフで行う限り不要だが、
  効率化のために使う場合の値）。
- 5.5 節の `<scf.SO.factor` に、入力に現れない原子種がある場合の扱い（既存実装は原子種の数だけ行を読む）。
- 部分内殻補正の密度の格子外での扱い（交換相関の章で決める）。
- 旧形式の実ファイルが手元にないため、旧形式の読み込みは記述どおりの合成ファイルでしか確認できていません。

## テスト

| 種類 | 内容 | 許容誤差の目安 |
|---|---|---|
| 解析 | s 型 Gauss 射影子 \\( e^{-\gamma r^2} \\) と s 型 Gauss 基底：\\( V\_{ab} = \varepsilon\\, S(\alpha,\gamma;\boldsymbol\tau\_{IK})\\, S(\gamma,\beta;\boldsymbol\tau\_{KJ}) \\)（二中心積分の章の閉じた式）。三次元格子での直接積分とも比較 | \\( 10^{-10} \\)（格子比較は \\( 10^{-6} \\)） |
| 解析 | p・d 型 Gauss 射影子と基底：Gauss 積の定理による重なりの積 | \\( 10^{-10} \\) |
| 解析 | MBK の \\( B^{-1} \\) 形式と対角化した \\( \pm 1 \\) 形式の一致 | 丸め誤差 |
| 解析 | \\( L\_x, L\_y, L\_z \\)（実基底）が \\( \mathbf L^2 = l(l+1) \\)、\\( [L\_x, L\_y] = iL\_z \\) を満たす | \\( 10^{-12} \\) |
| 解析 | \\( P\_\pm \\) が冪等・相補的で、トレースが \\( 2j+1 \\)、Clebsch–Gordan スピノルの射影と一致（\\( l = 1,2,3 \\)） | \\( 10^{-10} \\) |
| 解析 | 減衰関数の値と 2 階微分までの連続性 | 厳密 |
| 不変量 | Hermite 性：スカラーで \\( V\_{ab}(\mathbf L) = V\_{ba}(-\mathbf L) \\)、スピン軌道ありで 4 ブロック全体が Hermite | 丸め誤差 |
| 不変量 | 時間反転：\\( (i\sigma\_y)\\, \overline{V}\\, (i\sigma\_y)^{\dagger} = V \\)（スピン軌道ありの行列全体） | 丸め誤差 |
| 不変量 | スピン軌道ありで \\( \varepsilon\_{p+} = \varepsilon\_{p-} \\)、\\( b\_{p+} = b\_{p-} \\) のとき、スカラーの結果と一致し \\( V^{\uparrow\downarrow} = 0 \\) | 丸め誤差 |
| 不変量 | \\( \tfrac12 \mathrm{Tr}\_\sigma \\) が 4 節の演算子平均に一致 | 丸め誤差 |
| 不変量 | 構造全体の回転：スカラーの行列は実球面調和関数の回転行列で、スピン軌道ありの行列はさらにスピンの SU(2) 回転で変換される | \\( 10^{-10} \\) |
| 不変量 | 第三の原子 \\( K \\) を並進しても他の原子と一緒に動かせば不変、力の総和が 0 | \\( 10^{-10} \\) |
| 不変量 | 解析微分と有限差分の一致（第三の原子の微分、減衰関数の区間を含む） | \\( 10^{-6} \\) 程度 |
| 照合 | VPS 読み込み：`C_PBE19.vps`（j 依存、PCC あり）、`Pt_Sch.vps`（スカラー）、`C02.vps`（エネルギーが \\( \pm1 \\) 以外）、`C5Si5.vps`（最小ヘッダー）、`E.vps`（射影子 0 個）で、射影子数・\\( l \\)・エネルギー・列数・任意の \\( r \\) での値を fixture と比較 | 相対 \\( 10^{-10} \\) |
| 照合 | j 平均の結果（\\( \bar b\_p \\)、\\( \bar\varepsilon\_p \\)）を `Ba_PBE19.vps` のように符号が異なる射影子を含むファイルで確認 | 相対 \\( 10^{-12} \\) |
| 照合 | 小さな分子（例：CH₄、Pt₂）の非局所行列 \\( V^{\mathrm{NL}}(\mathbf R) \\) を既存実装の出力と比較（スカラー、スピン軌道あり、`scf.SO.factor` あり） | 要決定（\\( 10^{-8} \\) Hartree 程度を目安） |

検算スクリプト：`uv run --with numpy --with scipy --with sympy python tools/check_nonlocal.py`
