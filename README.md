# Granule Polarity Analysis (分泌顆粒極性解析プログラム)

このリポジトリは、HE染色組織画像における分泌顆粒の極性（Polarity Index）を自動解析するためのPythonスクリプトを公開しています。  
HE染色の生物学的特性（核＝ヘマトキシリン紫、顆粒＝エオジンピンク）に基づき、核シグナルを除外してエオジン（好酸性顆粒）特異的に輝度加重中心を算出する処理を実装しています。

---

## 概要 (Overview)

本プログラムは、HE染色組織画像から主要領域（細胞・腺腔等）を検出し、以下の2つの中心（Centroid）の位置ズレを基に**極性インデックス（Polarity Index）**を算出します。

1. **幾何学的中心 (Geometric Centroid):** 領域全体の形状の中心
2. **エオジン顆粒輝度重心 (Intensity-Weighted Centroid):** HSV変換とGreenチャンネル反転による、エオジン特異的輝度重心（核シグナル除外済み）

$$ \text{Polarity Index} = \frac{\text{重心間距離 (Distance)}}{\text{等価円相当径 (Equivalent Diameter)}} $$

---

## 主な機能 (Key Features)

* **エオジン特異的抽出:** HSV色空間を用いた好酸性シグナルの分離およびヘマトキシリン（核）シグナルの完全除外
* **自動領域抽出:** Otsuの大津二値化を用いた組織領域の自動セグメンテーション
* **ノイズ除去:** 一定面積（デフォルト: 5000 pixels）以下の微小ノイズの自動除外
* **定量的解析:** 面積、重心間距離、相当径、極性インデックスの自動算出
* **結果出力:** 解析結果の画面表示および CSV ファイル（`polarity_analysis_results.csv`）への自動保存

---

## 動作環境・必要ライブラリ (Requirements)

* **Python:** 3.8 以上推奨

### 必須ライブラリ
* `opencv-python` (cv2)
* `numpy`
* `pandas`
* `tabulate`

```bash
pip install opencv-python numpy pandas tabulate
