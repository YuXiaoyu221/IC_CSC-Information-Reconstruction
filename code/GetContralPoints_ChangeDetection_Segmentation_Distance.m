function [Updated_I,Updated_Imask]=GetContralPoints_ChangeDetection_Segmentation_Distance(Img_I0,I0_Mask,Img_I1,I1_Mask,SLIC_CloudArea,Kmeans_All,str_path,maxm,iwindow_s,thresholds)
%%% Fill some missing pixels before information reconstruction to control the coloring effects caused by subsequent isophote fusion. 
%%% A color-structure consistency constraint is constructed to solve these pixels with high reliability.
%%% When more than one pixel satisfies the color-structure consistency constraint, select the pixels cloest to the current missing pixel
%%% target image: Img_I0; Reference image: Img_I1
iwindow_s=floor(iwindow_s/2);% the window size for comparing local similarities: iwindow_s*2+1
thresholds_abs=thresholds(1);% the threshold to determine the similarity between I_r(p) and I_r(q) (1<t<3)
thresholds_cor=thresholds(2);% the threshold to determine the consistency between I_r(p) and I_t(p) (0<t<1)
if (min(maxm,size(SLIC_CloudArea,3))==0) 
    Updated_I=Img_I0;
    Updated_Imask=I0_Mask;
    return; 
end
Img_I0=double(Img_I0);inputData1=double(Img_I1);I0_Mask=double(I0_Mask);I1_Mask=double(I1_Mask);
inputData0_1=Img_I0;DataMask0_1=I0_Mask;
[m,n,bs]=size(Img_I0);
ChangeArea=I0_Mask+I1_Mask;
% maxm=size(SLIC_CloudArea,3);
%% definition of some parameters
slide_window=2;% enlarge the search window by 2
nKnownPoints=0;
iwindow_SS=iwindow_s*2+1;%the local window = iwindow_s*2+1
threshold_Sp=bs;%% when sum(abs(p-q))<=threshold_Sp,p and q are similiar pixels
threshold_Wp_ssd=(iwindow_s*2+1)*(iwindow_s*2+1)*bs*thresholds_abs;% when sum(abs(N_p-N_q))<=threshold_Wp_ssd,N_p and N_p are similiar
threshold_Wp_cor=thresholds_cor;%%when abs(corr(N_p_I0)-corr(N_p_I1))<=threshold_Wp_cor,N_p_I0 and N_p_I1 are consistency
%% creat dataset for similiar pixels searching
%%%Reference Image:J(m,n,b)--> J1(m,n,b*iwindow_SS*iwindow_SS)
%%%J1(i,j,:)is the neighborhoods of J(i,j) with in iwindow_SS*iwindow_SS
% inputData0_windows=zeros(m,n,bs*iwindow_SS*iwindow_SS);%target image
inputData1_windows=zeros(m,n,bs*iwindow_SS*iwindow_SS);%reference image
InDataMask0_Windows=255+zeros(m,n,1);%% if the nerghbor pixels of (i,j) are cloud pixels, it is not available
InDataMask0_Windows(iwindow_s+1:m-iwindow_s,iwindow_s+1:n-iwindow_s,:)=0;
for i=1:iwindow_SS %%left-up:£¨1,1£©start
    for j=1:iwindow_SS %%right-bottom:£¨iwindow_SS£¬iwindow_SS£©end
        for ibb=1:bs
        temp_data=inputData1(i:i+m-iwindow_SS,j:j+n-iwindow_SS,ibb);% temp_data=zeros(m+1-iwindow_SS,n+1-iwindow_SS,bs);
        inputData1_windows(iwindow_s+1:m-iwindow_s,iwindow_s+1:n-iwindow_s,(ibb-1)*iwindow_SS*iwindow_SS+(j-1)*iwindow_SS+i)=temp_data;
%         temp_data=inputData0(i:i+m-iwindow_SS,j:j+n-iwindow_SS,ibb);% temp_data=zeros(m+1-iwindow_SS,n+1-iwindow_SS,bs);
%         inputData0_windows(iwindow_s+1:m-iwindow_s,iwindow_s+1:n-iwindow_s,(ibb-1)*iwindow_SS*iwindow_SS+(j-1)*iwindow_SS+i)=temp_data;
        end
        temp_data=I0_Mask(i:i+m-iwindow_SS,j:j+n-iwindow_SS,:);
        InDataMask0_Windows(iwindow_s+1:m-iwindow_s,iwindow_s+1:n-iwindow_s,:)=temp_data+ InDataMask0_Windows(iwindow_s+1:m-iwindow_s,iwindow_s+1:n-iwindow_s,:);
    end
end
inputData1_windows_sum=sum(inputData1_windows,3);
inputData1_windows_1=zeros(m*n,1,bs*iwindow_SS*iwindow_SS);
for ibb=1:bs*iwindow_SS*iwindow_SS
    aa=inputData1_windows(:,:,ibb);
    inputData1_windows_1(:,:,ibb)=aa(:);
end
InDataMask0_Windows(InDataMask0_Windows>0)=1;
InDataMask0_Windows=1-InDataMask0_Windows;%% 0:cloud area; cloud-free area:1
ChangeArea(ChangeArea>0)=1;%% 
ChangeArea=1-ChangeArea;%% un-available area:0£¬available area:1
InDataMask0_Windows=InDataMask0_Windows.*ChangeArea;%% 
clear temp_data inputData1_windows
%% index matrix
ind_matrix=1:m*n;
ind_matrix=reshape(ind_matrix,m,n,1);
%% the correlation coefficient between target image and reference image
Corr_Images=zeros(m,n);
cw=iwindow_s+1;cw=3;
if (exist([str_path,'\Corr_Images_',num2str(cw),'.mat'])==0)
for i=cw+1:m-cw %row
    for j=cw+1:n-cw %coloum
        if(InDataMask0_Windows(i,j)==0)  continue;   end
        X0=Img_I0(i-cw:i+cw,j-cw:j+cw,:); X0=X0(:);
        Y0=inputData1(i-cw:i+cw,j-cw:j+cw,:); Y0=Y0(:);
        Corr_mnt=sum(X0.*Y0)/(sqrt(sum(X0.*X0))*sqrt(sum(Y0.*Y0))+0.00001);
        Corr_Images(i,j)=Corr_mnt;
    end
end
save([str_path,'\Corr_Images_',num2str(cw),'.mat'],'Corr_Images');
else
    Corr_Images=load([str_path,'\Corr_Images_',num2str(cw),'.mat']);
    Corr_Images=Corr_Images.Corr_Images;
end
%% search similiar pixels using segmatation results========
nkk=size(Kmeans_All,3)-1;% the total classes of Kmeans segmatation result
SLIC_CloudArea([1:iwindow_s,m-iwindow_s:m],:,:)=0;%%Prevents out-of-image
SLIC_CloudArea(:,[1:iwindow_s,n-iwindow_s:n],:)=0;%%Prevents out-of-image
for nci=1:min(maxm,size(SLIC_CloudArea,3))%% the counts of cloud areas in image
    I_Mask1=(SLIC_CloudArea(:,:,nci));%%
    ibl_n=max(max(I_Mask1));% the total numbers of superpixels
    [icloud_h,icloud_w]=find(I_Mask1>0);
    if (isempty(icloud_h)) continue; end
    nci_hs=min(icloud_h)-1;nci_he=max(icloud_h)+1;nci_ws=min(icloud_w)-1;nci_we=max(icloud_w)+1;
    nci_hs=max(nci_hs,0);nci_he=min(nci_he,m);nci_ws=max(nci_ws,0);nci_we=min(nci_we,n);
    nci_m=nci_he-nci_hs+1; nci_n=nci_we-nci_ws+1;
    ind_matrix_nci=ind_matrix(nci_hs:nci_he,nci_ws:nci_we,:);
    Kmeans_All_nci=double(Kmeans_All(nci_hs:nci_he,nci_ws:nci_we,:));
    inputData1_nci=inputData1(nci_hs:nci_he,nci_ws:nci_we,:);
    inputData1_windows_sum_nci=inputData1_windows_sum(nci_hs:nci_he,nci_ws:nci_we,:);
    inputData1_windows_1_nci=inputData1_windows_1(ind_matrix_nci,:,:);
    ind_matrix_nci=1:nci_m*nci_n;
    I_Mask1_nci=double(I_Mask1(nci_hs:nci_he,nci_ws:nci_we));
    DataMask0_1_nci=DataMask0_1(nci_hs:nci_he,nci_ws:nci_we);inputData0_1_nci=inputData0_1(nci_hs:nci_he,nci_ws:nci_we,:);
    nKnownPoints_t=0;
    InDataMask0_SearchBufferArea_0=imfilter(double(I_Mask1),ones(iwindow_s,iwindow_s)/iwindow_s/iwindow_s);
    InDataMask0_SearchBufferArea_0(InDataMask0_SearchBufferArea_0>0)=1;    
    %%%%% ===get the cloest one
    InDataMask0_SearchBufferArea_t0=imfilter(InDataMask0_SearchBufferArea_0,ones(slide_window*2,slide_window*2)/slide_window/slide_window/4);
    InDataMask0_SearchBufferArea_t0(InDataMask0_SearchBufferArea_t0>0)=1;
    InDataMask0_SearchPoints=InDataMask0_SearchBufferArea_t0-InDataMask0_SearchBufferArea_0;
    SearchP_H=ind_matrix(InDataMask0_SearchPoints==1);%%
    [SearchP_H,SearchP_W]=ind2sub([m,n],SearchP_H);
    while isempty(SearchP_H)==0 && nKnownPoints_t<ibl_n %%
       for ipp=1:length(SearchP_H) 
           if Corr_Images(SearchP_H(ipp),SearchP_W(ipp))<threshold_Wp_cor continue; end
           ipp_class=Kmeans_All(SearchP_H(ipp),SearchP_W(ipp),nkk+1);%%the class of current pixel
           ClassMask_ipp=Kmeans_All_nci(:,:,ipp_class);%1: the same class; 0: different class
           ClassMask_ipp=ClassMask_ipp.*(I_Mask1_nci);%%%%
           if (any(ClassMask_ipp(:)>0)==0) continue; end
           temRef_pp_ij=inputData1(SearchP_H(ipp),SearchP_W(ipp),:);%%
           temRef_pp_ij_mat=repmat(temRef_pp_ij,nci_m,nci_n,1);%%            
           dif_ij=sum(abs(inputData1_nci-temRef_pp_ij_mat),3);
           dif_ij(ClassMask_ipp==0)=threshold_Sp+threshold_Sp;
           if (isempty(find(dif_ij<threshold_Sp, 1))) continue; end
           %%get the neighboring pixels£ºtemRef_pp_windows
           temRef_pp_windows=(inputData1(SearchP_H(ipp)-iwindow_s:SearchP_H(ipp)+iwindow_s,SearchP_W(ipp)-iwindow_s:SearchP_W(ipp)+iwindow_s,:));
           temRef_pp_windows_lines(1,1,:)=temRef_pp_windows(:);
           inputData1_windows_sum_dif=abs(inputData1_windows_sum_nci-sum(temRef_pp_windows_lines));
           inputData1_windows_sum_dif(dif_ij>=threshold_Sp)=threshold_Wp_ssd+threshold_Wp_ssd;
           %%%==when sum(abs(N_p-N_q))<=threshold_Wp_ssd,N_p and N_p are similiar
           ind_DifHW=ind_matrix_nci(inputData1_windows_sum_dif<=threshold_Wp_ssd);                  
           if isempty(ind_DifHW)  continue;  end                
           inputData1_windows_t=inputData1_windows_1_nci(ind_DifHW,:,:);%%% when sum(abs(p-q))<=threshold_Sp,p and q are similiar pixels
           temRef_pp_windows_lines_mat=repmat(temRef_pp_windows_lines,length(ind_DifHW),1,1);
           inputData1_windows_sum_dif=sum(abs(inputData1_windows_t-(temRef_pp_windows_lines_mat)),3);    
           if isempty(find(inputData1_windows_sum_dif<=threshold_Wp_ssd, 1)) continue; end           
           ind_DifHW=ind_DifHW(inputData1_windows_sum_dif<=threshold_Wp_ssd);%%
           inputData1_windows_sum_dif=inputData1_windows_sum_dif(inputData1_windows_sum_dif<=threshold_Wp_ssd);
           [ind_DifH,ind_DifW]=ind2sub([nci_m,nci_n],ind_DifHW);
           %%%%When more than one pixel£¬ select the pixels cloest to the current missing pixel
           m_kpoints_lable=I_Mask1_nci(ind_DifHW);
           for ikp=1:length(ind_DifH)
               %%%% updated I_Mask1
               ibl_i=I_Mask1_nci(ind_DifH(ikp),ind_DifW(ikp));
               if (ibl_i==0)  continue; end
               I_Mask1_nci(I_Mask1_nci==ibl_i)=0;          
               SamePoints_ind=m_kpoints_lable(m_kpoints_lable==ibl_i);            
               if (length(SamePoints_ind)>1) %% select the nearest one
                   inputData1_windows_sum_dif_s=inputData1_windows_sum_dif(m_kpoints_lable==ibl_i);
%                    dif_ij_s=dif_ij(ind_DifHW(SamePoints_ind));
%                    min_dif_windows_s=dif_ij_s+inputData1_windows_sum_dif_s/(iwindow_s*2+1)/(iwindow_s*2+1);                 
                   [min_dif_indh,min_dif_indw]=ind2sub([nci_m,nci_n],ind_DifHW(inputData1_windows_sum_dif_s==min(inputData1_windows_sum_dif_s)));                   
                   %%% updated image I
                   DataMask0_1_nci(min_dif_indh(1),min_dif_indw(1))=0;
                   inputData0_1_nci(min_dif_indh(1),min_dif_indw(1),:)=Img_I0(SearchP_H(ipp),SearchP_W(ipp),:);
               elseif (length(SamePoints_ind)==1)
                   %%%updated image I
                   DataMask0_1_nci(ind_DifH(ikp),ind_DifW(ikp))=0;
                   inputData0_1_nci(ind_DifH(ikp),ind_DifW(ikp),:)=Img_I0(SearchP_H(ipp),SearchP_W(ipp),:);
               end
               nKnownPoints_t=nKnownPoints_t+1;
           end                 
       end
       %%%%% ===until all superpixels are finished
       InDataMask0_SearchBufferArea_tt=imfilter(InDataMask0_SearchBufferArea_t0,ones(slide_window*2,slide_window*2)/slide_window/slide_window/4);
       InDataMask0_SearchBufferArea_tt(InDataMask0_SearchBufferArea_tt>0)=1;
       InDataMask0_SearchPoints=InDataMask0_SearchBufferArea_tt-InDataMask0_SearchBufferArea_t0;
       SearchP_H=ind_matrix(InDataMask0_SearchPoints==1);
       [SearchP_H,SearchP_W]=ind2sub([m,n],SearchP_H);
       InDataMask0_SearchBufferArea_t0=InDataMask0_SearchBufferArea_tt;
    end
    %%% updating
    DataMask0_1_t=zeros(m,n);DataMask0_1_t(nci_hs:nci_he,nci_ws:nci_we,:)=DataMask0_1_nci;
    I_Mask1(I_Mask1>0)=1; DataMask0_1_t(DataMask0_1_t>0)=1;
    DataMask0_1_t=double(I_Mask1)-DataMask0_1_t;
    DataMask0_1(DataMask0_1_t>0)=0;
    inputData0_1_t=zeros(m,n,bs);inputData0_1_t(nci_hs:nci_he,nci_ws:nci_we,:)=inputData0_1_nci;
    DataMask0_1_t=repmat(DataMask0_1_t,1,1,bs);
    inputData0_1(DataMask0_1_t>0)=inputData0_1_t(DataMask0_1_t>0);
    nKnownPoints=nKnownPoints+nKnownPoints_t;
end
 Updated_I=uint8(inputData0_1);
 Updated_Imask=uint8(DataMask0_1);
end