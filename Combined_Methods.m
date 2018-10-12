clear all;
close all;
clc;

%Input image
img = imread ('01_dr.JPG');
img = rgb2gray(img);
img = im2double (img);

%% Canny
%Value for Thresholding
T_Low = 0.075;
T_High = 0.175;

%Gaussian Filter Coefficient
B = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5;4, 9, 12, 9, 4;2, 4, 5, 4, 2 ];
B = 1/159.* B;

%Convolution of image by Gaussian Coefficient
A=conv2(img, B, 'same');

%Filter for horizontal and vertical direction
KGx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
KGy = [1, 2, 1; 0, 0, 0; -1, -2, -1];

%Convolution by image by horizontal and vertical filter
Filtered_X = conv2(A, KGx, 'same');
Filtered_Y = conv2(A, KGy, 'same');

%Calculate directions/orientations
arah = atan2 (Filtered_Y, Filtered_X);
arah = arah*180/pi;

pan=size(A,1);
leb=size(A,2);

%Adjustment for negative directions, making all directions positive
for i=1:pan
    for j=1:leb
        if (arah(i,j)<0) 
            arah(i,j)=360+arah(i,j);
        end;
    end;
end;

arah2=zeros(pan, leb);

%Adjusting directions to nearest 0, 45, 90, or 135 degree
for i = 1  : pan
    for j = 1 : leb
        if ((arah(i, j) >= 0 ) && (arah(i, j) < 22.5) || (arah(i, j) >= 157.5) && (arah(i, j) < 202.5) || (arah(i, j) >= 337.5) && (arah(i, j) <= 360))
            arah2(i, j) = 0;
        elseif ((arah(i, j) >= 22.5) && (arah(i, j) < 67.5) || (arah(i, j) >= 202.5) && (arah(i, j) < 247.5))
            arah2(i, j) = 45;
        elseif ((arah(i, j) >= 67.5 && arah(i, j) < 112.5) || (arah(i, j) >= 247.5 && arah(i, j) < 292.5))
            arah2(i, j) = 90;
        elseif ((arah(i, j) >= 112.5 && arah(i, j) < 157.5) || (arah(i, j) >= 292.5 && arah(i, j) < 337.5))
            arah2(i, j) = 135;
        end;
    end;
end;

%figure, imagesc(arah2); colorbar;

%Calculate magnitude
magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
magnitude2 = sqrt(magnitude);

BW = zeros (pan, leb);

%Non-Maximum Supression
for i=2:pan-1
    for j=2:leb-1
        if (arah2(i,j)==0)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i,j+1), magnitude2(i,j-1)]));
        elseif (arah2(i,j)==45)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j-1), magnitude2(i-1,j+1)]));
        elseif (arah2(i,j)==90)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j), magnitude2(i-1,j)]));
        elseif (arah2(i,j)==135)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j+1), magnitude2(i-1,j-1)]));
        end;
    end;
end;

BW = BW.*magnitude2;
%figure, imshow(BW);

%Hysteresis Thresholding
T_Low = T_Low * max(max(BW));
T_High = T_High * max(max(BW));

T_res = zeros (pan, leb);

for i = 1  : pan
    for j = 1 : leb
        if (BW(i, j) < T_Low)
            T_res(i, j) = 0;
        elseif (BW(i, j) > T_High)
            T_res(i, j) = 1;
        %Using 8-connected components
        elseif ( BW(i+1,j)>T_High || BW(i-1,j)>T_High || BW(i,j+1)>T_High || BW(i,j-1)>T_High || BW(i-1, j-1)>T_High || BW(i-1, j+1)>T_High || BW(i+1, j+1)>T_High || BW(i+1, j-1)>T_High)
            T_res(i,j) = 1;
        end;
    end;
end;

edge_final = uint8(T_res.*255);
%Show final edge detection result
% figure, imshow(edge_final);
% title('Canny Edge Detection How?');
figure, imagesc(edge_final);colormap(gray); 
title ('Canny Edge Detection How?');

%% local variance based edge detection
I = img;
% define the window size
sz = 3;
mn = floor(sz/2);

% preallocate the matrix
output = zeros(size(I));

% pad the matrix with zeros
I = padarray(I,[mn mn]);

for i = 1:size(I,1)-mn*2
    for j = 1:size(I,2)-mn*2
        tmp = I(i:i+(sz-1),j:j+(sz-1));
        mu = mean(tmp(:));
        tmp2 = mean(tmp(:).^2);
        output(i,j) = tmp2 - mu.^2;
    end
end
figure, imagesc(uint8(output));colormap(gray); 
title ('Local Varience Based Edge Detection');

%% Robert
b = img;
[m,n]=size(img);
L(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        L(i,j)=-1*b(i,j)+0+0+1*b(i+1,j+1);
    end;
end;



M(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        M(i,j)=0-1*b(i,j+1)+1*b(i+1,j)+0;
    end;
end;

%%figure;
%subplot(2,2,1)
%imshow(L)
%title('Robert    Gx');
%subplot(2,2,2)
%imshow(M)
%title('Robert    Gy');
N=M+L;
%subplot(2,2,3)
%imshow(N)
%title('Robert    Gx+Gy');
%subplot(2,2,4)
%imshow(b)
%%title('Original Image');
figure, imagesc(uint8(N));colormap(gray); 
title ('Robert');

%% PREWIT
N(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        N(i,j)=-1*b(i,j)-1*b(i,j+1)-1*b(i,j+2)+0+0+0+1*b(i+2,j)+1*b(i+2,j+1)+1*b(i+2,j+2);
    end;
end;
O(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        O(i,j)=-1*b(i,j)+0+1*b(i,j+2)-1*b(i+2,j)+0+1*b(i+1,j+2)-1*b(i+2,j)+0+1*b(i+2,j+2);
    end;
end;

%figure;
%subplot(2,2,1)
%imshow(N)
%title('Prewit    Gx');
%subplot(2,2,2)
%imshow(O)
%title('Prewit    Gy');
Z=N+O;
%subplot(2,2,3)
%imshow(Z)
%title('Prewit    Gx+Gy');
%subplot(2,2,4)
%imshow(b)
%title('Original Image');

figure, imagesc(uint8(Z));colormap(gray); 
title ('Prewit');

%% SOBEL
P(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        P(i,j)=-1*b(i,j)-2*b(i,j+1)-1*b(i,j+2)+0+0+0+1*b(i+2,j)+2*b(i+2,j+1)+1*b(i+2,j+2);
    end;
end;


R(1:m,1:n)=0;
for i=1:m-2;
    for j=1:n-2;
        R(i,j)=-1*b(i,j)+0+1*b(i,j+2)-2*b(i+1,j)+0+2*b(i+1,j+2)-1*b(i+2,j)+0+1*b(i+2,j+2);
    end;
end;

 %figure;
%subplot(2,2,1)
%imshow(P)
%title('Sobel    Gx');
%subplot(2,2,2)
%imshow(R)
%title('Sobel   Gy');

Y=P+R;
%subplot(2,2,3)
%imshow(Y)
%title('Soble  Gx+Gy');
%subplot(2,2,4)
%imshow(b)

%title('Original  Image');

figure, imagesc(uint8(Y));colormap(gray); 
title ('Sobel');

%% Laplacian and Gaussian (LoG)
hsize = [7 7]; %only runs for odd integer step_size
sigma = 0.3;
find_edges(img,hsize,sigma) 

%% kirsch Edge Detection
img_k = kirschedge(img);

%% Robinson Edge Detection
img_r = robinsonedge(img);


%% Marr hildreth edge detector

%smoothening the image with a filter
gfilter= [0 0 1 0 0;
       0 1 2 1 0;
       1 2 -16 2 1;
       0 1 2 1 0;
       0 0 1 0 0];
   
smim=conv2(img,gfilter);

% finding the zero crossings

[rr,cc]=size(smim);
zc=zeros([rr,cc]);

for i=2:rr-1
    for j=2:cc-1
        if (smim(i,j)>0)
             if (smim(i,j+1)>=0 && smim(i,j-1)<0) || (smim(i,j+1)<0 && smim(i,j-1)>=0)
                             
                zc(i,j)= smim(i,j+1);
                        
            elseif (smim(i+1,j)>=0 && smim(i-1,j)<0) || (smim(i+1,j)<0 && smim(i-1,j)>=0)
                    zc(i,j)= smim(i,j+1);
            elseif (smim(i+1,j+1)>=0 && smim(i-1,j-1)<0) || (smim(i+1,j+1)<0 && smim(i-1,j-1)>=0)
                  zc(i,j)= smim(i,j+1);
            elseif (smim(i-1,j+1)>=0 && smim(i+1,j-1)<0) || (smim(i-1,j+1)<0 && smim(i+1,j-1)>=0)
                  zc(i,j)=smim(i,j+1);
            end
                        
        end
            
    end
end


otpt=im2uint8(zc);
% thresholding
otptth= otpt>105;

% figure;
%   subplot(2,2,1);imshow(im);title('Origional image');
%   subplot(2,2,2);imshow(smim);title('Smoothened image');
%   subplot(2,2,3);imshow(otpt);title('Output image');
%  subplot(2,2,4);imshow(otptth);title('Output image with threshold');

  % final result
figure, imagesc(uint8(otptth));colormap(gray); 
title ('Marr Hildreth edge detector AKA LoG');
   


