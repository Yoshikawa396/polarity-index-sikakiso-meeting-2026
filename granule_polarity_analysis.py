"""
Granule Polarity Analysis Script
Analyzes the polarity of secretory granules in tissue images.
"""

import cv2
import numpy as np
import pandas as pd
import os

def analyze_granule_polarity(image_path: str, label: str):
    """
    Analyzes the polarity of secretory granules in a tissue image.
    Calculates the Polarity Index based on the distance between 
    the Intensity Weighted Centroid (Eosinophilic Granules) and the Geometric Centroid (Cell/Region).
    """
    # 1. Load Image in BGR Color Space (修正: グレースケールではなくカラーで読み込み)
    img_bgr = cv2.imread(image_path)
    if img_bgr is None:
        raise FileNotFoundError(f"Image not found at path: {image_path}")

    # -------------------------------------------------------------------------
    # 【修正点1】エオジン（好酸性顆粒）特異的シグナルの抽出
    # HE染色ではエオジン（ピンク）はGreen(G)チャンネルを強く吸収するため、
    # 255 - Green チャンネルをとることで顆粒シグナルを強調。
    # さらにHSV変換を用いてヘマトキシリン（核＝青紫）の領域を除外するマスクを作成。
    # -------------------------------------------------------------------------
    img_hsv = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2HSV)
    
    # 好酸性（ピンク〜赤系）のHSV範囲定義 (Hue: 140-175付近)
    lower_eosin = np.array([140, 30, 50])
    upper_eosin = np.array([175, 255, 255])
    eosin_mask = cv2.inRange(img_hsv, lower_eosin, upper_eosin)

    # Greenチャンネルを反転して好酸性輝度マップを作成
    img_g = img_bgr[:, :, 1]
    eosin_intensity = 255 - img_g
    
    # エオジンマスクを適用し、核や背景のシグナルを完全カット（顆粒のみの輝度情報）
    granule_signal = cv2.bitwise_and(eosin_intensity, eosin_intensity, mask=eosin_mask)

    # -------------------------------------------------------------------------
    # 【修正点2】組織/細胞領域（関心領域: ROI）のセグメンテーション
    # -------------------------------------------------------------------------
    img_gray = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2GRAY)
    _, thresh = cv2.threshold(img_gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    # 輪郭抽出
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    results = []

    for i, cnt in enumerate(contours):
        area = cv2.contourArea(cnt)
        
        # 間質ノイズや微小領域を除外 (閾値設定)
        if area < 5000:
            continue

        # 現在の対象領域のバイナリマスク作成
        region_mask = np.zeros_like(img_gray)
        cv2.drawContours(region_mask, [cnt], -1, 255, -1)

        # 2. Calculate Geometric Centroid (組織/細胞領域の幾何学的中心)
        M_geom = cv2.moments(cnt)
        if M_geom["m00"] == 0:
            continue
        cx_geom = M_geom["m10"] / M_geom["m00"]
        cy_geom = M_geom["m01"] / M_geom["m00"]

        # 3. Calculate Intensity-Weighted Centroid (【修正点3】エオジン顆粒特異的輝度重心)
        # 抽出した顆粒シグナルに対し、ROIマスクを適用
        granule_masked = cv2.bitwise_and(granule_signal, granule_signal, mask=region_mask)
        M_intensity = cv2.moments(granule_masked)
        
        if M_intensity["m00"] == 0:
            continue
        cx_weighted = M_intensity["m10"] / M_intensity["m00"]
        cy_weighted = M_intensity["m01"] / M_intensity["m00"]

        # 4. Calculate Distance (幾何学的中心と顆粒輝度重心のシフト距離)
        distance = np.sqrt((cx_weighted - cx_geom)**2 + (cy_weighted - cy_geom)**2)

        # 5. Calculate Equivalent Diameter (相当直径)
        equiv_diameter = np.sqrt(4 * area / np.pi)

        # 6. Calculate Polarity Index (規格化された極性インデックス)
        polarity_index = distance / equiv_diameter if equiv_diameter > 0 else 0

        # 結果の保存
        results.append({
            "Image_ID": f"HE_{label}",
            "Cell_or_Acinus_ID": f"Region_{len(results)+1:02d}",
            "Area_pixels": round(area, 1),
            "Distance_pixels": round(distance, 2),
            "Equiv_Diameter": round(equiv_diameter, 2),
            "Polarity_Index": round(polarity_index, 3),
            "Time_Point": label
        })

    return pd.DataFrame(results)


# ==========================================
# 実行例
# ==========================================
if __name__ == "__main__":
    image_files = {
        "0 min": "photo_0min.jpg",
        "200 min": "photo_200min.jpg",
        "300 min": "photo_300min.jpg"
    }

    all_results = []

    for time_label, path in image_files.items():
        if os.path.exists(path):
            print(f"Analyzing {time_label}...")
            df_res = analyze_granule_polarity(path, time_label)
            all_results.append(df_res)
        else:
            print(f"Warning: File '{path}' not found. Skipping...")

    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        print("\n=== Analysis Completed ===")
        print(final_df.to_markdown())
        
        final_df.to_csv("polarity_analysis_results.csv", index=False)
        print("Results saved to 'polarity_analysis_results.csv'.")
