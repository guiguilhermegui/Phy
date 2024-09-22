function segmented_image = edge_segmentation(image)

    % Convert to double precision for numerical operations
    image = double(image);

    % Detect edges using the Canny edge detector
    edges = edge(image, 'canny');

    % Fill in the interiors of detected regions
    segmented_image = imfill(edges, 'holes');

end
