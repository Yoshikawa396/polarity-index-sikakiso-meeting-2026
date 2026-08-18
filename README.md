# Granule Polarity Analysis (分泌顆粒極性解析プログラム)

このリポジトリは、組織画像における分泌顆粒の極性（Polarity Index）を自動解析するためのPythonスクリプトを公開しています。  
本プログラムは学会発表（2026年）に関連する解析コードです。

---

## 概要 (Overview)

本プログラムは、HE染色等の組織画像から主要領域（細胞・腺腔等）を検出し、以下の2つの中心（Centroid）の位置ズレを基に**極性インデックス（Polarity Index）**を算出します。

1. **幾何学的重心 (Geometric Centroid):** 領域全体の形状の中心
2. **輝度加重重心 (Intensity-Weighted Centroid):** 顆粒や核などの濃密領域（輝度）を考慮した重心

$$ \text{Polarity Index} = \frac{\text{重心間距離 (Distance)}}{\text{等価円相当径 (Equivalent Diameter)}} $$

---

## 主な機能 (Key Features)

* **自動領域抽出:** Otsuの大津二値化（Otsu's thresholding）を用いた組織領域の自動セグメンテーション
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
pip install opencv-python numpy pandas tabulat``` polarity-index-sikakiso-meeting-2026 polarity-index-sikakiso-meeting-2026 polarity-index-sikakiso-meeting-2026 polarity-index-sikakiso-meeting-2026 polarity-index-sikakiso-meeting-2026 polarity-index-sikakiso-meeting-2026 polarity-index-sikakiso-meeting-2026
