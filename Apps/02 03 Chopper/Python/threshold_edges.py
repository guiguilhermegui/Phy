import cv2
import numpy as np

# Load the image
img = cv2.imread('image.jpg', 0)

# Apply edge detection and thresholding
edges = cv2.Canny(img, 100, 200)
_, thresh = cv2.threshold(edges, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

# Find connected components
connectivity = 8  # 8-connected components
num_labels, labels = cv2.connectedComponents(thresh, connectivity)

# Create subfolders
output_dir = 'shapes'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Save individual shapes
for label in range(1, num_labels):
    mask = np.zeros_like(img)
    mask[labels == label] = 255
    shape_img = cv2.bitwise_and(img, mask)
    shape_name = f'shape_{label}.jpg'
    cv2.imwrite(os.path.join(output_dir, shape_name), shape_img)
