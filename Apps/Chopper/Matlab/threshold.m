function segmented_image = threshold_segmentation(image)

    % Convert to double precision for numerical operations
    image = double(image);

    % Determine the threshold using Otsu's method
    level = graythresh(image);

    % Apply thresholding
    segmented_image = im2bw(image, level);

end
