clc
clear
close all

%%
[file, path] = uigetfile({'*.jpg;*.png'},'Please select an image:');
FileName     = fullfile(path, file);
Im = imread(FileName);

%%
Im_Gray = double(im2gray(Im));
[m, n]  = size(Im_Gray);
Im_Edge_X = abs(diff(Im_Gray,1,1));
Im_Edge_Y = abs(diff(Im_Gray,1,2));
Im_Edge_Cnt_X = zeros(1,m-1);
Im_Edge_Cnt_Y = zeros(1,n-1);
for ii = 1:m-1
    idx_temp = find(Im_Edge_X(ii,:)>0);
    Im_Edge_Cnt_X(ii) = sum(diff(idx_temp(1:end-1))==1 & diff(idx_temp,2)==0);
end
for ii = 1:n-1
    idx_temp = find(Im_Edge_Y(:,ii)>0);
    Im_Edge_Cnt_Y(ii) = sum(diff(idx_temp(1:end-1))==1 & diff(idx_temp,2)==0);
end
Max_X = max(Im_Edge_Cnt_X) - median(Im_Edge_Cnt_X);
Max_Y = max(Im_Edge_Cnt_Y) - median(Im_Edge_Cnt_Y);
[peak_x_edge, idx_x_edge] = findpeaks(Im_Edge_Cnt_X,"MinPeakProminence",Max_X/2);
[peak_y_edge, idx_y_edge] = findpeaks(Im_Edge_Cnt_Y,"MinPeakProminence",Max_Y/2);
temp_left  = mean(Im_Edge_Cnt_X(idx_x_edge+(-5:-3).'),1);
temp_right = mean(Im_Edge_Cnt_X(idx_x_edge+(3:5).'),1);
temp_idx   = peak_x_edge < max(temp_left,temp_right)+Max_X/4;
idx_x_edge(temp_idx)  = [];
peak_x_edge(temp_idx) = [];
temp_left  = mean(Im_Edge_Cnt_Y(idx_y_edge+(-5:-3).'),1);
temp_right = mean(Im_Edge_Cnt_Y(idx_y_edge+(3:5).'),1);
temp_idx   = peak_y_edge < max(temp_left,temp_right)+Max_Y/4;
idx_y_edge(temp_idx)  = [];
peak_y_edge(temp_idx) = [];

%%
figure
subplot(5,5,[1:4 6:9 11:14 16:19])
imagesc(Im)
if size(Im,3)==1
    colormap gray
end
axis off
subplot(5,5,21:24)
plot(1:n-1,Im_Edge_Cnt_Y)
hold on
plot(idx_y_edge,peak_y_edge,'r^')
axis([-inf inf -inf inf])
box off
set(gca,'XTick',[],'YTick',[])
subplot(5,5,5:5:20)
plot(Im_Edge_Cnt_X,1:m-1)
hold on
plot(peak_x_edge,idx_x_edge,'r>')
axis([-inf inf -inf inf])
box off
set(gca,'ydir','reverse','XTick',[],'YTick',[])
drawnow

%%
idx_x_edge = [0 idx_x_edge m];
idx_y_edge = [0 idx_y_edge n];
X_Num = length(idx_x_edge)-1;
Y_Num = length(idx_y_edge)-1;
Border_Val = zeros(length(idx_x_edge)-1,length(idx_y_edge)-1,4);
for ii = 1:X_Num
    xx = idx_x_edge(ii)+1:idx_x_edge(ii+1);
    for jj = 1:Y_Num
        yy = idx_y_edge(jj)+1:idx_y_edge(jj+1);
        temp_data = Im_Gray(idx_x_edge(ii)+2,yy);
        Border_Val(ii,jj,1) = mean(temp_data(:));
        temp_data = Im_Gray(xx,idx_y_edge(jj+1)-2);
        Border_Val(ii,jj,2) = mean(temp_data(:));
        temp_data = Im_Gray(idx_x_edge(ii+1)-2,yy);
        Border_Val(ii,jj,3) = mean(temp_data(:));
        temp_data = Im_Gray(xx,idx_y_edge(jj)+2);
        Border_Val(ii,jj,4) = mean(temp_data(:));
    end
end
Border_Val(isnan(Border_Val)) = 0;
H = zeros(X_Num*Y_Num,X_Num*Y_Num);
for ii = 1:X_Num
    for jj = 1:Y_Num
        if ii<X_Num
            idx1 = sub2ind([X_Num,Y_Num],ii,jj);
            idx2 = sub2ind([X_Num,Y_Num],ii+1,jj);
            H(idx1,idx1) = H(idx1,idx1) + Border_Val(ii,jj,3)^2;
            H(idx2,idx2) = H(idx2,idx2) + Border_Val(ii+1,jj,1)^2;
            H(idx1,idx2) = H(idx1,idx2) - 2*Border_Val(ii,jj,3)*Border_Val(ii+1,jj,1);
        end
        if jj<Y_Num
            idx1 = sub2ind([X_Num,Y_Num],ii,jj);
            idx2 = sub2ind([X_Num,Y_Num],ii,jj+1);
            H(idx1,idx1) = H(idx1,idx1) + Border_Val(ii,jj,2)^2;
            H(idx2,idx2) = H(idx2,idx2) + Border_Val(ii,jj+1,4)^2;
            H(idx1,idx2) = H(idx1,idx2) - 2*Border_Val(ii,jj,2)*Border_Val(ii,jj+1,4);
        end
    end
end
Lambda = 10;
H = (H+H.')/2;
W = (0.01+(0:X_Num*Y_Num-1)).';
H = H + Lambda*eye(size(H)).*W;
A_ans = quadprog(H,-2*Lambda*W,[],[],[],[],zeros(X_Num*Y_Num,1));
Block_Scale_Total = reshape(A_ans,X_Num,Y_Num);

%%
Im_Adjust = zeros(size(Im));
for ii = 1:length(idx_x_edge)-1
    xx = idx_x_edge(ii)+1:idx_x_edge(ii+1);
    for jj = 1:length(idx_y_edge)-1
        yy = idx_y_edge(jj)+1:idx_y_edge(jj+1);
        Im_Adjust(xx,yy,:) = Block_Scale_Total(ii,jj) * double(Im(xx,yy,:));
    end
end

%%
if size(Im_Adjust,3)>1
    R = Im_Adjust(:,:,1);
    G = Im_Adjust(:,:,2);
    B = Im_Adjust(:,:,3);
    Im_Adjust_Gray = 0.2989 * R + 0.5870 * G + 0.1140 * B;
else
    Im_Adjust_Gray = Im_Adjust;
end
Min_Im_Gray = quantile(Im_Gray(:),0.001);
Max_Im_Gray = quantile(Im_Gray(:),0.999);
Min_Im_Adjust_Gray = quantile(Im_Adjust_Gray(:),0.001);
Max_Im_Adjust_Gray = quantile(Im_Adjust_Gray(:),0.999);
Im_Adjust = (Im_Adjust - Min_Im_Adjust_Gray) / (Max_Im_Adjust_Gray - Min_Im_Adjust_Gray);
Im_Adjust = Im_Adjust * (Max_Im_Gray - Min_Im_Gray) + Min_Im_Gray;
Im_Adjust = uint8(Im_Adjust);

%%
figure
set(gcf,'WindowState', 'maximized')
subplot(1,2,1)
imshow(Im)
title('Original Iamge')
subplot(1,2,2)
h = imshow(Im);
title('Image blocks adjustment...')
drawnow
for ii = 1:length(idx_x_edge)-1
    xx = idx_x_edge(ii)+1:idx_x_edge(ii+1);
    for jj = 1:length(idx_y_edge)-1
        yy = idx_y_edge(jj)+1:idx_y_edge(jj+1);
        h.CData(xx,yy,:) = Im_Adjust(xx,yy,:);
        drawnow
    end
end
title('Done!')

%%
imwrite(Im_Adjust,[path '\' file(1:end-4) '-Adjusted' file(end-3:end)])