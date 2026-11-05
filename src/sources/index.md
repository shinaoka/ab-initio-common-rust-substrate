# 公開資料の一覧

実装者が直接読んでよい資料です。URL は 2026-11-05 時点で確認しました。

## マニュアル

- OpenMX Ver. 4.0 マニュアル：[HTML](https://www.openmx-square.org/openmx_man4.0/index.html)、
  [PDF](https://www.openmx-square.org/openmx_man4.0/openmx4.0.pdf)、
  [日本語版](https://www.openmx-square.org/openmx_man4.0jp/index.html)
- 再実装に特に関係する節：
  - 11 汎関数、12 基底関数（12.1 概要、12.2 primitive basis、12.4 PAO の最適化）、13 擬ポテンシャル
  - 14 カットオフエネルギーと格子、15 k 点、16 SCF の収束
  - 27 O(N) 法、40 DFT+U、45 交換相互作用、46 NEGF、47 Wannier 関数
  - 70 後処理コード向けのインターフェース（70.1 scfout ファイル）、71 OpenMX API
- ADPACK Ver. 2.2 マニュアル：[HTML](https://www.openmx-square.org/adpack_man2.2/adpack2_2.html)

## 技術ノート

索引：<https://www.openmx-square.org/tech_notes/tech_notes.html>

| 題目 | URL | このノートの関連章 |
|---|---|---|
| Total Energy and Forces (Ver. 1.3) | [tech1-1_3](https://www.openmx-square.org/tech_notes/tech1-1_3/tech1-1_3.html) | 二中心積分、全エネルギーと力 |
| Stress tensor | [tech23](https://www.openmx-square.org/tech_notes/tech23-1_0.pdf) | 全エネルギーと力 |
| Efficient Projector Expansion in LCAO | [tech15](https://www.openmx-square.org/tech_notes/tech15-1_0/tech15-1_0.html) | 中性原子ポテンシャル |
| Numerical Atomic Basis Orbitals from H to Kr | [tech14](https://www.openmx-square.org/tech_notes/tech14-1_0/tech14-1_0.html) | PAO 基底 |
| Non-Collinear Spin Density Functional | [tech2](https://www.openmx-square.org/tech_notes/tech2-1_0/tech2-1_0.html) | （未作成） |
| O(N) Krylov Subspace Method | [tech3](https://www.openmx-square.org/tech_notes/tech3-1_1/tech3-1_1.html) | （未作成） |
| O(N) Divide-Conquer with LNO | [tech26](https://www.openmx-square.org/tech_notes/tech26-1_0/tech26-1_0.html) | （未作成） |
| Efficient Integration of Green's Functions | [tech4](https://www.openmx-square.org/tech_notes/tech4-1_1/tech4-1_1.html) | 密度行列 |
| Orbital Optimization Method | [tech5](https://www.openmx-square.org/tech_notes/tech5-1_0/tech5-1_0.html) | PAO 基底 |
| LDA+U Method | [tech6](https://www.openmx-square.org/tech_notes/tech6-1_0/tech6-1_0.html) | （未作成） |
| Closest Wannier Function | [cwf](https://www.openmx-square.org/tech_notes/cwf/cwf.html) | （未作成） |
| Geometry Optimization | [tech7](https://www.openmx-square.org/tech_notes/tech7-1_0/tech7-1_0.html) | （未作成） |
| Extrapolation of Charge Density | [tech8](https://www.openmx-square.org/tech_notes/tech8-1_0/tech8-1_0.html) | 密度行列と電荷ミキシング |
| Charge Analysis | [tech9](https://www.openmx-square.org/tech_notes/tech9-1_0/tech9-1_0.html) | （未作成） |
| Householder Method for Tridiagonalization | [tech10](https://www.openmx-square.org/tech_notes/tech10-1_0/tech10-1_0.html) | 一般化固有値 |
| Electric Polarization by Berry Phase | [tech11](https://www.openmx-square.org/tech_notes/tech11-1_1/tech11-1_1.html) | （未作成） |
| Charge Mixing Methods | [tech12](https://www.openmx-square.org/tech_notes/tech12-1_0/tech12-1_0.html) | 密度行列と電荷ミキシング |
| Two Electron Repulsion Integrals | [exx/JCP2009](https://www.openmx-square.org/exx/JCP2009/index.html) | （未作成） |
| Fast Spherical Bessel Transform | [exx/CPC2009](https://www.openmx-square.org/exx/CPC2009/index.html) | 二中心積分 |
| Electronic Transport Calculations | [tech13](https://www.openmx-square.org/tech_notes/tech13-1_0/tech13-1_0.html) | （未作成） |
| Eigenchannels and Current density | [tech20](https://www.openmx-square.org/tech_notes/tech20-1_0/tech20-1_0.html) | （未作成） |
| Low-order Scaling Method | [tech16](https://www.openmx-square.org/tech_notes/tech16-1_0/tech16-1_0.html) | （未作成） |
| Finite Element Method for Atoms | [tech17](https://www.openmx-square.org/tech_notes/tech17-1_0/tech17-1_0.html) | （未作成） |
| Unfolding Optical Transition Weights | [PDF](https://www.openmx-square.org/tech_notes/unfolding_otw.pdf) | （未作成） |
| Absolute Binding Energies of Core Levels | [tech24](https://www.openmx-square.org/tech_notes/tech24-1_0/tech24-1_0.html) | （未作成） |
| Calculation of Optical Matrix Elements | [PDF](https://www.openmx-square.org/tech_notes/Opt_MatElements.pdf) | （未作成） |
| Dielectric Function and Conductivity | [PDF](https://www.openmx-square.org/tech_notes/Dielectric_Function_YTL.pdf) | （未作成） |
| Spin texture in k-space | [PDF](https://www.openmx-square.org/tech_notes/note_kSpin-1_0.pdf) | （未作成） |
| Nudged Elastic Band Method | [PDF](https://www.openmx-square.org/tech_notes/NEB.pdf) | （未作成） |

## 論文

マニュアル 3 節「Related papers」より。

| 内容 | 文献 |
|---|---|
| PAO 基底の変分最適化 | T. Ozaki, PRB **67**, 155108 (2003) |
| H〜Kr の数値原子基底 | T. Ozaki and H. Kino, PRB **69**, 195113 (2004) |
| 射影子展開 | T. Ozaki and H. Kino, PRB **72**, 045121 (2005) |
| 高速球 Bessel 変換 | J. Toyoda and T. Ozaki, CPC **181**, 277 (2010) |
| Krylov 部分空間 O(N) 法 | T. Ozaki, PRB **74**, 245101 (2006) |
| 低次スケーリング法 | T. Ozaki, PRB **82**, 075131 (2010) |
| DC-LNO 法 | T. Ozaki, M. Fukuda, G. Jiang, PRB **98**, 245137 (2018) |
| 大規模並列化 | T. V. T. Duy and T. Ozaki, CPC **185**, 777 (2014); CPC **185**, 153 (2014) |
| DFT+U | M. J. Han, T. Ozaki, J. Yu, PRB **74**, 045110 (2006) |
| 交換相互作用 | M. J. Han, T. Ozaki, J. Yu, PRB **70**, 184421 (2004) |
| NEGF | T. Ozaki, K. Nishio, H. Kino, PRB **81**, 035116 (2010) |
| ESM 法 | T. Ohwaki et al., JCP **136**, 134101 (2012) |
| MLWF | H. Weng, T. Ozaki, K. Terakura, PRB **79**, 235118 (2009) |
| CWF | T. Ozaki, PRB **110**, 125115 (2024) |
| バンドアンフォールディング | C.-C. Lee, Y. Yamada-Takamura, T. Ozaki, JPCM **25**, 345501 (2013) |
| 内殻準位の束縛エネルギー | T. Ozaki and C.-C. Lee, PRL **118**, 026401 (2017) |
| 精度比較（Δ ゲージ） | K. Lejaeghere et al., Science **351**, aad3000 (2016) |

一般的な手法（Kleinman–Bylander 形式、Blöchl の多重射影子、Troullier–Martins 擬ポテンシャル、
Morrison–Bylander–Kleinman 形式、Pulay/Kerker ミキシングなど）は原論文を参照します。
