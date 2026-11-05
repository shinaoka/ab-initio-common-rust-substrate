# H(R)・S(R) 疎ブロック形式〔設計〕

> この章は設計の叩き台です。OpenMX 内部のデータ配置を写すものではなく、
> 新しい実装のインターフェースとして定めます。

## 対象

局在基底の実空間行列

$$
O\_{I\alpha, J\beta}(\mathbf L) = \langle \phi\_{I\alpha}(\mathbf r - \mathbf R\_I) | \hat O | \phi\_{J\beta}(\mathbf r - \mathbf R\_J - \mathbf L) \rangle
$$

を、原子対と格子ベクトルごとの密ブロックとして保持します。\\( \hat O \\) は重なり（恒等演算子）、
運動エネルギー、Hamiltonian などです。

## 構造と値の分離

- **構造**：原子対 \\( (I, J, \mathbf L) \\) のリスト、各原子の軌道数、ブロックの形状、値の配列内の位置。
  原子配置・基底・カットオフが変わらない限り不変です。
- **値**：ブロックの数値。SCF の反復ごとに更新されます。

構造を一度作り、値だけを差し替えられるようにします。

## 決めるべき規約

| 項目 | 叩き台 |
|---|---|
| 原子対の範囲 | \\( \lVert \mathbf R\_J + \mathbf L - \mathbf R\_I \rVert \le r\_{c,I} + r\_{c,J} \\) を満たすすべての組（自分自身を含む。既存実装と同じく境界を含む。境界上の対の値は 0） |
| 対の重複 | \\( (J, I, -\mathbf L) \\) も保持する（Hermite 性で復元できるが、利用側の簡単さを優先）。ただし O(N) 解法（`DC`、`Krylov`、`DC-LNO`）の密度行列は上の意味で Hermite でなく \\( \rho\_{I\alpha,J\beta}(\mathbf L) \ne \rho\_{J\beta,I\alpha}(-\mathbf L) \\) なので、復元はできず**両方の向きを必ず保持する**（[O(N) 法](../modules/on-methods.md)） |
| 格子ベクトルの表現 | 整数の 3 つ組 \\( (n\_1, n\_2, n\_3) \\) |
| ブロック内の並び | 行は \\( I \\) の軌道、列は \\( J \\) の軌道。軌道の並びは [PAO 基底](../modules/pao-basis.md) の規約 |
| メモリー配置 | 未定（列優先を第一候補。faer・tenferro との受け渡しを考慮） |
| スピン | 共線：スピンごとに別の行列。非共線：2×2 のスピンブロックを持つ複素行列 |
| Bloch 和 | \\( O(\mathbf k) = \sum\_{\mathbf L} e^{+i\mathbf k\cdot\mathbf L} O(\mathbf L) \\)。位相は格子ベクトルだけに付け、原子位置は含めない（既存実装と同じ。下の〔規約〕） |
| k 点の表現 | 逆格子ベクトルに対する分数座標 \\( \mathbf q \\)。\\( \mathbf k\cdot\mathbf L = 2\pi\\,\mathbf q\cdot\mathbf n \\) |
| 非共線の k 空間行列 | \\( 2N\times 2N \\)、スピンが外側・軌道が内側（\\( p = \sigma N + a \\)）。重なりは \\( I\_2\otimes S(\mathbf k) \\) |

## 既存実装の k 空間の規約〔規約〕

設計の選択を既存実装と揃えるため、確認した事実を記録します（詳細は [k 点と Bloch 和](../modules/kpoints-bloch.md)）。

- 位相は \\( e^{+2\pi i\\,\mathbf q\cdot\mathbf n} \\) で、\\( \mathbf n \\) は列側（ket 側）の原子 \\( J \\) の周期像を指す格子ベクトルの整数座標です。
  原子位置は位相に入りません（技術ノート 1 の Eq. 1、28–29 とも一致）。同じ \\( (I, J) \\) の複数の周期像は加算されます。
- Cluster（Γ 点）では \\( \sum\_{\mathbf L} O(\mathbf L) \\) を実数の行列として作ります。
- 非共線では \\( \uparrow\uparrow \\)、\\( \downarrow\downarrow \\)、\\( \uparrow\downarrow \\) の 3 ブロックを Bloch 和で作り、
  \\( \downarrow\uparrow \\) ブロックは \\( \uparrow\downarrow \\) ブロックの Hermite 共役で埋めます。
  実空間では、各スピンブロックの実部と、スピン軌道相互作用などによる虚部を別々に保持しています。
- 実空間の密度行列は \\( \rho(\mathbf L) = \sum\_{\mathbf k} w\_{\mathbf k}\sum\_\mu f\_{\mu\mathbf k}\\, e^{+i\mathbf k\cdot\mathbf L}\\, \overline{c\_{I\alpha,\mu}}\\, c\_{J\beta,\mu} \\)
  （共線では実部）で、\\( H(\mathbf L) \\) と同じ位相の向きです。

## 未決事項〔未確定〕

- 非共線の場合の実部・虚部、スピン成分の保持方法（既存実装は実部と虚部を別の配列に分けて持つが、新しい形式では複素数のブロックにするか）。
- 共線計算で、時間反転による k 点の削減を利用側（密度行列の組み立て）に任せるか、形式の側で \\( \mathrm{Re} \\) をとる関数を提供するか。
- 原子位置を含む位相の規約（後処理で必要になる場合）を、変換関数として提供するかどうか。
- MPI 分散時に、どの原子対をどのランクが持つか。
