"""
Granule Polarity Analysis Script
Analyzes the polarity of secretory granules in tissue images.
"""

import os
import cv2
import numpy as np
import pandas as pd


def analyze_granule_polarity(image_path: str, label: str):
    """
    Analyzes the polarity of secretory granules in a tissue image.
    Calculates the Polarity Index based on the distance between 
    the Intensity Weighted Centroid (Granules) and the Geometric Centroid (Cell/Region).
    """
    # 1. Load Image in Grayscale
    img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise FileNotFoundError(f"Image not found at path: {image_path}")

    # 2. Invert image intensity so that darker regions (granules/nuclei) have higher values
    inv_img = 255 - img

    # 3. Apply Otsu's thresholding to segment the tissue/cell region
    _, thresh = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    # 4. Find contours of the main tissue/cell regions
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    results = []

    for i, cnt in enumerate(contours):
        area = cv2.contourArea(cnt)
        
        # Filter out small noisy regions (Threshold can be adjusted)
        if area < 5000:
            continue

        # Create a binary mask for the current region
        mask = np.zeros_like(img)
        cv2.drawContours(mask, [cnt], -1, 255, -1)

        # 5. Calculate Geometric Centroid (Unweighted Center of the Region)
        M_geom = cv2.moments(cnt)
        if M_geom["m00"] == 0:
            continue
        cx_geom = M_geom["m10"] / M_geom["m00"]
        cy_geom = M_geom["m01"] / M_geom["m00"]

        # 6. Calculate Intensity-Weighted Centroid (Granule Center of Mass)
        intensity_masked = cv2.bitwise_and(inv_img, inv_img, mask=mask)
        M_intensity = cv2.moments(intensity_masked)
        
        if M_intensity["m00"] == 0:
            continue
        cx_weighted = M_intensity["m10"] / M_intensity["m00"]
        cy_weighted = M_intensity["m01"] / M_intensity["m00"]

        # 7. Calculate Distance (Shift between Geometric and Weighted Centroid)
        distance = np.sqrt((cx_weighted - cx_geom)**2 + (cy_weighted - cy_geom)**2)

        # 8. Calculate Equivalent Diameter
        equiv_diameter = np.sqrt(4 * area / np.pi)

        # 9. Calculate Polarity Index (Normalized Distance)
        polarity_index = distance / equiv_diameter if equiv_diameter > 0 else 0

        # Store results
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
# Execution Example
# ==========================================
if __name__ == "__main__":
    # Define image file paths and corresponding time points
    # Replace these paths with your local image files
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

    # Combine results if images were processed
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        print("\n=== Analysis Completed ===")
        print(final_df.to_markdown())
        
        # Save to CSV
        final_df.to_csv("polarity_analysis_results.csv", index=False)
        print("Results saved to 'polarity_analysis_results.csv'.")