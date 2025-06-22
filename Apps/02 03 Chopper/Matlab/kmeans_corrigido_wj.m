clc
clear

%readfolder = '/MATLAB Drive/RDI';
%readfolder = 'G:\0DD\Desktop\Pastas\VSCODE\unb\Phy\Apps\02 03 Chopper\Matlab';
readfolder = 'G:\0DD\Desktop\Workspace\Phy\partisan\pics\RDI_ArcelorMittal_SEM\0_RDI';
file = '9.51-2-EDS.tif';
splitted = split(file, '.');
filename = string(splitted(2));
filetype = string(splitted(3));
folder = filename;
medfiltsz = [10 10];
numColors = 3;
NumAttempts = 1;

cd(readfolder)
A = imread(file);

%% Equalização de Histograma
H = histeq(A);

%% Conversão em tons de cinza
G = im2gray(H);

%% Filtro Mediana
F = medfilt2(G, medfiltsz);

%% Gradiente (para watershed)
grad = imgradient(F);

%% Segmentação K-means
L = imsegkmeans(A, numColors);
B = labeloverlay(A, L);
figure();
imshowpair(L, B, 'montage'); title("Labeled Image Overlay RGB");
%% === Watershed Controlado por Marcadores ===

% 1. Marcadores dos objetos (foreground markers)
se = strel('disk', 20);
opened = imopen(F, se);
fgm = imregionalmax(opened);

% Opcional: Limpar pequenos objetos dos marcadores
fgm = bwareaopen(fgm, 20);

% 2. Marcadores do fundo (background markers)
bw = imbinarize(F);
D = bwdist(~bw);
DL = watershed(D);
bgm = DL == 0;

% 3. Impor mínimos no gradiente
grad2 = imimposemin(grad, bgm | fgm);

% 4. Watershed controlado por marcadores
W = watershed(grad2);

% Visualização dos marcadores (opcional)
figure;
imshowpair(fgm, bgm, 'montage');
title('Marcadores: Foreground (esq) e Background (dir)');

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
clusterDisplayAndWatershed(A, B, ab, numColors, NumAttempts)


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

function clusterDisplayAndWatershed(A, B, ab, numColors, NumAttempts)

% Inputs:
%   A: Original image (RGB or grayscale).
%   B: Image to extract clusters from (e.g., grayscale or a filtered version of A).
%   ab: a*b* color space representation of the image.
%   numColors: Number of clusters (k).
%   NumAttempts: Number of attempts for k-means.

% Perform k-means clustering
pixel_labels = imsegkmeans(ab,numColors,NumAttempts=NumAttempts);

% Overlay labels on the original image
B2 = labeloverlay(A, pixel_labels);
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


function cropAndSaveImages(A, stats, filename)
for i = 1:qtd
    bbox = stats(i).BoundingBox;  % Extract the bounding box

    % Check if the bounding box is valid (within image bounds)
    bbox(1) = max(1, bbox(1));
    bbox(2) = max(1, bbox(2));
    bbox(3) = min(size(A,2) - bbox(1) +1, bbox(3));
    bbox(4) = min(size(A,1) - bbox(2) +1, bbox(4));

    if bbox(3) <= 0 || bbox(4) <= 0
        warning('Bounding box %d is invalid and will be skipped.', i);
        continue; % Skip to the next iteration if bbox is invalid
    end

    CHOPPER = imcrop(A, bbox); % Crop the original image using the bounding box

    formatSpec = 'file_%s_%d.tif';
    imwrite(CHOPPER, sprintf(formatSpec,filename,i), 'BitDepth', 16);  % Save the cropped image

end

end