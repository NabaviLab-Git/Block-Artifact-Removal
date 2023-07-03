clc
clear
close all

%%
Im = imread('cameraman.tif');
figure
subplot(1,2,1)
imshow(Im)
m = size(Im,1);
n = size(Im,2);
M = 20;
N = 15;
idx_x = sort([0,round(4*randn(1,M)+linspace(m/M,m-m/M,M)),m],'ascend');
idx_y = sort([0,round(4*randn(1,N)+linspace(n/N,n-n/N,N)),n],'ascend');
s_factor = rand(M+1,N+1)*0.5+0.5;
for ii = 1:length(idx_x)-1
    xx = idx_x(ii)+1:idx_x(ii+1);
    for jj = 1:length(idx_y)-1
        yy = idx_y(jj)+1:idx_y(jj+1);
        Im(xx,yy,:) = uint8(double(Im(xx,yy,:))*s_factor(ii,jj));
    end
end
subplot(1,2,2)
imshow(Im)