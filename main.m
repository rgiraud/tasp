% This code is free to use for any non-commercial purposes.
% It contains an implementation of the SCALP superpixel method proposed in:
% [1] - Texture-Aware Superpixel Segmentation
%       Rémi Giraud, Vinh-Thong Ta, Nicolas Papadakis, Yannick Berthoumieu
%       IEEE International Conference on Image Processing (ICIP 2019)
%
% Note that the core of the implementation is based on the provided code associated to the following paper:
% [2] - Zhengqin Li, Jiansheng Chen
%       Superpixel Segmentation using Linear Spectral Clustering
%       International Conference on Computer Vision and Pattern Recognition (CVPR), 2015
%
% If you use this code, please cite both [1] and [2].
%
% (C) Rémi Giraud, 2019
% rgiraud@u-bordeaux.fr, https://remi-giraud.enseirb-matmeca.fr/
% Bordeaux-INP, IMS Laboratory

clear; close all; clc

addpath(genpath('utils_sp'))

%Compilation
mex -O CFLAGS="\$CFLAGS -Wall" ./TASP_mex.cpp -outdir ./

% Image loading
data_path = './data';
img_name = 'texture_compo.png'; 
% img_name = 'test_img.jpg';
I = imread(sprintf('%s/%s',data_path, img_name));
[h,w,z] = size(I);

SP_nbr = 90; 
compactness = 0.1;
S = TASP_mex(uint8(I), SP_nbr, compactness); 

% Plotting results
B = ones(h,w,3);
for i=1:max(S(:))
    j = (S == i);
    bb = bwboundaries(j);
    if ~isempty(bb)
        for k=1:length(bb{1})
            B(bb{1}(k,1),bb{1}(k,2),1:2) = 0;
            B(bb{1}(k,1),bb{1}(k,2),3) = 255;
        end
    end
end

figure,
imagesc(double(I)/255.*B); 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%To use display and evaluation superpixel functions
%%download the toolbox at https://remi-giraud.enseirb-matmeca.fr/research/#superpixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


