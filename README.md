# OpenMX 実装ノート

[OpenMX](https://www.openmx-square.org/) の手法を、OpenMX のソースコードを読まずに再実装できる水準まで記述する実装仕様書です。
[mdBook](https://rust-lang.github.io/mdBook/) で HTML にレンダリングできます。

## ビルド

```sh
cargo install mdbook   # 未導入の場合
mdbook serve --open
```

## 構成

- `src/process/` 運用規則（クリーンルームの手順、記述の規則、ライセンス）
- `src/sources/` 公開資料の索引と、公開資料の空白
- `src/design/` 記法と全体設計
- `src/modules/` モジュールごとの仕様
- `src/testing/` テストの方針
- `tools/` 式の検算用スクリプト

## ライセンス

- 文章：[CC BY 4.0](LICENSE-CC-BY-4.0.txt)
- コード例・スクリプト・テストデータ：[MIT](LICENSE-MIT)

OpenMX 本体と VPS・PAO データベース（GPL）のコード・データは含みません。
詳しくは `src/process/license.md` を参照してください。
