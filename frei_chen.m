%% Frei-Chen Approach for edge detection
clear all
close all
clc

%A = [55, 63, 77, 80, 129;22,94,87,27,26;107,154,76,30,37;108,252,220,75,30;20,111,185,51,10];
%A  = imread('cameraman.tif');
A = imread('Satellite-Image-Port-Cape-Town-South-Africa.jpg');
A = rgb2gray(A);
A = im2double(A);
k = [-1,0,-1;-sqrt(2),0,-sqrt(2);-1,0,1];

Hsl = convolve(A, k);
figure,
%Hsl = uint8(Hsl);
imshow(Hsl);
title ('First Approach');

%% Approach 2
% pad the initial matrix with zero  
[t s] = size (k);
B = padarray(A,[t-2 s-2]);

% PRE-ALLOCATE THE MATRIX
Output = zeros([size(A,1) size(A,2)]);

% PERFORM COONVOLUTION
for i = 1:size(B,1)-2;
    for j = 1:size(B,2)-2;
        Temp = B(i:i+2,j:j+2).*k;
        disp(Temp);
        Output(i,j) = sum(Temp(:));
    end
 end
figure,
%Output = uint8(Output);
imshow(Output);
title('Second Approach');

%% Approach 3
b = A;
[m,n]=size(A);
L(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        L(i,j)=-1*b(i,j)+0-1*b(i,j+2)-sqrt(2)*b(i+1,j)+0-sqrt(2)*b(i+1,j+2)-1*b(i+2,j)+0+1*b(i+2,j+2);
    end;
end;
figure, 
%L = uint8(L);
imshow(L);
title ('Third Approach');