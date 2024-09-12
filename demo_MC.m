clear all
warning('off')

% CBM3D
addpath('qtfm/') % 

rate = 75;

% gt
O_Img = imread('./starfish.png');

% mask
mask = im2double(imread('75_mask.png'));
        
% Sampling
[m,n,c] = size(O_Img);
D = mask .* double(O_Img);

% G  
figure
[G_NMF, ~] = inexact_alm_NMF_MC(D, mask, sqrt(m*n), eps, 1e-5, 600);
G_NMF = uint8(G_NMF);

imwrite(uint8(G_NMF), 'MC.png')


