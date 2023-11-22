function Kmeans_lable=ImgSegmatationRelable_Kmeans(KmeansData,nk)
%%% re-lable for KmeansData according to different classes
%%% Kmeans_lable: single band, 1¡Üp¡ÜKmeans_lable
[m,n,ibs]=size(KmeansData);
i_ssucess=0; nk_mat=zeros(1,ibs);
%%  re-lable for KmeansData according to different classes
Kmeans_lable0=zeros(m,n);
for ib=1:ibs
    %%%pixelValue: 0~255£¨x-coordinate£©; count of number£¨y-coordinate£©
    [count,~]=imhist(KmeansData(:,:,ib));
    indexi=find(count>0);
    if length(indexi)~=nk
        nk_mat(1,ib)=length(indexi);
        continue;
    end
    i_ssucess=1;
    Kmeans_lable0=double(KmeansData(:,:,ib));
    for ip=1:nk
        temp_p=indexi(ip)-1;%%
        Kmeans_lable0(Kmeans_lable0==temp_p)=-ip;
    end
end

%% special case
if (i_ssucess==0 && length(find(nk_mat==nk_mat(1)))==ibs)
    nk=nk_mat(1);
    for ib=1:ibs
        %%%pixelValue: 0~255£¨x-coordinate£©; count of number£¨y-coordinate£©
        [count,~]=imhist(KmeansData(:,:,ib));
        indexi=find(count>0);
        if length(indexi)~=nk
            nk_mat(1,ib)=length(indexi);
            continue;
        end
        Kmeans_lable0=double(KmeansData(:,:,ib));
        for ip=1:nk
            temp_p=indexi(ip)-1;
            Kmeans_lable0(Kmeans_lable0==temp_p)=-ip;
        end
    end    
end

%% write and output the result
Kmeans_lable0=uint8(-Kmeans_lable0);
Kmeans_lable=uint8(zeros(m,n,nk+1));
for ib=1:nk
    Kmeans_lablei=zeros(m,n);
    Kmeans_lablei(Kmeans_lable0==ib)=1;
    Kmeans_lable(:,:,ib)=Kmeans_lablei;
end
Kmeans_lable(:,:,nk+1)=Kmeans_lable0;
%%%  max(max(Kmeans_lable))
%%%  [count,pixelValue]=imhist(Kmeans_lable);
%%%  figure,imshow(Kmeans_lable(:,:,7)*20)
end