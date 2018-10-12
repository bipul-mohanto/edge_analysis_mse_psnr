clear all;
close all;
clc

I = imread('cameraman.tif');
x = [0 size(I,2)]; % horizontal direction 0 to number of column
% why this???
y = [size(I,1)/2 size(I,1)/2];
c = improfile(I,x,y);

figure
subplot(2,1,1)
imagesc(I);colormap('default'); 
hold on
plot(x,y,'r');
subplot(2,1,2);
%% for one channel
plot(c(:,1,1),'r');
hold on
%% for g and b channel
% plot(c(:,1,2),'g');
% plot(c(:,1,3),'b');
