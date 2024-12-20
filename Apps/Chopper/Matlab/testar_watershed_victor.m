% copiar e colar código https://www.mathworks.com/help/images/marker-controlled-watershed-segmentation.html
%testar até funfar
% unificar com o Chopper














function segmented_image = region_segmentation(image)

    % Convert to double precision for numerical operations
    image = double(image);

    % Apply morphological operations to enhance edges
    image = imtophat(image, strel('disk', 10));

    % Perform watershed segmentation
    [~, ~, segmented_image] = watershed(imcomplement(image));

end
