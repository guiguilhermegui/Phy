clc
clear

% Read the image
cd 'D:\Desktop\Pastas\Unb\partisan\Examples\RDI_ArcelorMittal_SEM\0102png'
filename = '9.51-2-EDS.png';

I = imread(filename);
% Convert the image to grayscale
%I_gray = rgb2gray(I);

% Apply a median filter to reduce noise
I_filtered = medfilt2(I);

% Calculate the gradient magnitude of the filtered image
gmag = imgradient(I_filtered);

% Determine an optimal threshold value using Otsu's method
level = graythresh(I_filtered);

% Binarize the image using the calculated threshold
bw = imbinarize(I_filtered, level);

% Fill holes in the binary image
bw = imfill(bw, 'holes');

% Apply the watershed transform to segment the image
L = watershed(gmag);

% Extract bounding boxes of segmented regions
stats = regionprops(L, 'BoundingBox');
cd ..
cd('PNG')
% Iterate through each region
for i = 1:numel(stats)
    % Extract the bounding box
    bbox = stats(i).BoundingBox;
    
    % Crop the original image using the bounding box
    cropped_img = imcrop(I, bbox);
    
    %discover the file name
    %[name,~,~] = fileparts(filename);

    % Save the cropped image
    %imwrite(cropped_img, sprintf('file_%s_%d.jpg',name,i));
    formatSpec ='file_%s_%d.jpg';
    imwrite(cropped_img, sprintf(formatSpec,filename,i));
end
