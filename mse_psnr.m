clear all;
close all;
clc;

img = imread ('Satellite-Image-Port-Cape-Town-South-Africa.jpg');
%img=  imread ('WorldView-3 Satellite Image Airport Mapping Madrid Spain.jpg');
img = rgb2gray(img);
img = im2double (img);
%%  Algo Section

cc = edge(img, 'sobel');
figure, imagesc(cc);colormap(gray); 
%title ('Local Varience Based Edge Detection');

%% MSE and PSNR

n = size(img);
M = n(1);
N = n(2);
MSE = sum(sum((img-cc).^2))/(M*N);

PSNR = 10*log10 (256*256/MSE);
fprintf('\n MSE: %7.5f', MSE);

fprintf('\n PSNR: %9.7f dB', PSNR);
