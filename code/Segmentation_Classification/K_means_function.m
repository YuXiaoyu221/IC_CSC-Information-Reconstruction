function [Img_segmentation,nk_color]=K_means_function(InputData,nk)
%%% Img_segmentation result
%%% nk_color 
[dim_i,dim_j,~] = size(InputData); % 
means = zeros(nk, 3); % Initialize means to randomly-selected colors in the original photo.
rand_x = ceil(dim_i*rand(nk, 1));%
rand_y = ceil(dim_j*rand(nk, 1));
for i = 1:nk
    means(i,:) = InputData(rand_x(i), rand_y(i), :);%the center of the class
end
for itr=1:100
    nk_color=zeros(nk,3);
    s_ind=zeros(nk,1);
    for i=1:dim_i
        for j=1:dim_j
            r=InputData(i,j,1);g=InputData(i,j,2);b=InputData(i,j,3);
            [~, ind]=min(sum((repmat([r,g,b],nk,1)-means).^2,2));
            %repmat(A,k,1)
            nk_color(ind,:)=nk_color(ind,:)+[r,g,b];
            s_ind(ind)=s_ind(ind)+1;
        end
    end
    for ii=1:nk
        if s_ind(ii)>0
            nk_color(ii,:)=nk_color(ii,:)./s_ind(ii);
        end
    end
    d=sum(sqrt(sum((nk_color-means).^2,2)));%distance
    if d<1e-5
        break
    end
    means=nk_color;   
end
means = round(means);
% itr
%% the color classes
% figure; hold on
% for i=1:k
%    col = (1/255).*means(i,:);
%    rectangle('Position', [i, 0, 1, 1], 'FaceColor', col, 'EdgeColor', col);
% end
% axis off
Img_segmentation =InputData; % double(imread('bird_large.tiff'));
for i = 1:dim_i
    for j = 1:dim_j
        r = Img_segmentation(i,j,1); g = Img_segmentation(i,j,2); b = Img_segmentation(i,j,3);
        [~, ind]=min(sum((repmat([r,g,b],nk,1)-means).^2,2));
        Img_segmentation(i,j,:) = means(ind,:);
    end 
end