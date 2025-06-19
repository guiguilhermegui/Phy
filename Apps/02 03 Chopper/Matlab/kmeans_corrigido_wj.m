clc
clear

readfolder = '/MATLAB Drive/RDI';
file = '9.51-2-EDS.tif';
splitted = split(file, '.');
filename = string(splitted(2)); 
filetype = string(splitted(3));
folder = filename; 
medfiltsz = [10 10]; 
numColors = 3; 
NumAttempts = 1;   

cd(readfolder)
O = imread(file); 

%% Equalização de Histograma
H = histeq(O); 

%% Conversão em tons de cinza
G = im2gray(H); 

%% Filtro Mediana
F = medfilt2(G, medfiltsz); 

%% Gradiente (para watershed)
grad = imgradient(F); 

%% Segmentação K-means
L = imsegkmeans(O, numColors); 
B = labeloverlay(O, L); 
figure();
imshowpair(L, B, 'montage'); title("Labeled Image Overlay RGB");

%% Watershed baseado no gradiente
W = watershed(grad); 


H_uint8 = im2uint8(H);
W_uint8 = im2uint8(W);
imshowpair(H_uint8, W_uint8, 'montage');
title("Original and Watershed");


%% Aplicando Kmeans Cluster

lab_B = rgb2lab(B);

cmap = parula(numColors); 
[~, ~, indexedImage] = unique(lab_B); 
ColorImage = ind2rgb(indexedImage, cmap); 

ab = lab_B(:,:,2:3);
ab = im2single(ab);

figure()
clusterDisplayAndWatershed(O, B, ab, numColors, NumAttempts)


%%  Aplicando a separacao por caixas

stats = regionprops(W, 'BoundingBox');

qtd = length(stats);
fprintf('%d partículas foram identificadas \n Aqui vai um preview das caixas de segmentação de partículas \n',qtd);


imshowpair(H,W,'montage');



   
    %% Exportando Imagens
    
    mkdir(folder);
    cd(folder);
    fprintf('exportando imagens para a pasta \n %s \n',folder)


fprintf('export terminado \n')

%% Funcoes 

function clusterDisplayAndWatershed(O, B, ab, numColors, NumAttempts)

% Inputs:
%   O: Original image (RGB or grayscale).
%   B: Image to extract clusters from (e.g., grayscale or a filtered version of O).
%   ab: a*b* color space representation of the image.
%   numColors: Number of clusters (k).
%   NumAttempts: Number of attempts for k-means.

% Perform k-means clustering
pixel_labels = imsegkmeans(ab,numColors,NumAttempts=NumAttempts);

% Overlay labels on the original image
B2 = labeloverlay(O, pixel_labels);
imshow(B2);
title("Labeled Image a*b*");

% Display each cluster
for i = 1:numColors
    mask = pixel_labels == i;
    cluster = B.*uint8(mask);
    figure; % Create a new figure for each cluster
    imshow(cluster);
    title(sprintf("Objects in Cluster %d", i));
end





end


function cropAndSaveImages(O, stats, filename)
for i = 1:qtd
    bbox = stats(i).BoundingBox;  % Extract the bounding box

    % Check if the bounding box is valid (within image bounds)
    bbox(1) = max(1, bbox(1));
    bbox(2) = max(1, bbox(2));
    bbox(3) = min(size(O,2) - bbox(1) +1, bbox(3));
    bbox(4) = min(size(O,1) - bbox(2) +1, bbox(4));

    if bbox(3) <= 0 || bbox(4) <= 0
        warning('Bounding box %d is invalid and will be skipped.', i);
        continue; % Skip to the next iteration if bbox is invalid
    end

    CHOPPER = imcrop(O, bbox); % Crop the original image using the bounding box

    formatSpec = 'file_%s_%d.tif';
    imwrite(CHOPPER, sprintf(formatSpec,filename,i), 'BitDepth', 16);  % Save the cropped image

end

end