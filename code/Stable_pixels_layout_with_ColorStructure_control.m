function [Updated_I,Updated_Imask,Updated_J,Updated_Jmask]=Stable_pixels_layout_with_ColorStructure_control(Img_I,I_mask,Img_J,J_mask,str_cell,p_StablePixelsSolving,p_SLICSeg,p_KmeansSeg)
%%% == Fill some missing pixels before information reconstruction to control the coloring effects caused by subsequent isophote fusion. 
%%% == A color-structure consistency constraint is constructed to solve these pixels with high reliability.
str_tif='.tif';
iflag_StablePixels=p_StablePixelsSolving{1,2}; % if "Img_I & I_mask" have been updated, it can be set as 1
iwindow_s=p_StablePixelsSolving{2,2};% the window size for comparing local similarities: iwindow_s¡Áiwindow_s£¨eg.3¡Á3£©
iflag_SLICSeg=p_SLICSeg{1,2}; % if iflag_Kmeans_SILC_Reable has been solved already, it can be set as 1
i_Superpixels=p_SLICSeg{2,2}; % the initial size of eack superpixel in SLIC segmatation result
iflag_KmeansSeg=p_KmeansSeg{1,2}; % if iflag_Kmeans_SILC_Reable has been solved already, it can be set as 1
i_KmeansClasses=p_KmeansSeg{2,2}; % the classes of Kmeans segmatation result
iflag_Kmeans_SILC_Reable=0; % if iflag_Kmeans_SILC_Reable has been solved already, it can be set as 1
str_InputPath=char(str_cell(1,3));
for TarImg=1:2
    if TarImg==1
        I=Img_I;  str_It=char(str_cell(1,1));  It_mask=I_mask;
        J=Img_J;  str_Jt=char(str_cell(1,2));  Jt_mask=J_mask;
    elseif TarImg==2
        I=Img_J;  str_It=char(str_cell(1,2));  It_mask=J_mask;
        J=Img_I;  str_Jt=char(str_cell(1,1));  Jt_mask=I_mask;
    end
    %%% the number of cloud areas in image
    L_Mask=bwlabel(uint8(It_mask),4);
    icpend=max(max(L_Mask));%%% the number of cloud areas in image
    if icpend>0
        %%%Stable_pixels_layout
%         iFlag_KnownPointPro=1; %% 1 means the result can be input directly
        if iflag_StablePixels==0
        % J_resample=imfilter(double(J),ones(3,3)/9);
        J_resample=J;
        Str_SLIC=[str_InputPath,'\',str_It,'_SLIC.png'];
        Str_Kmeans=[str_InputPath,'\',str_Jt,'_Kmeans.tif'];
%         iflag_Kmeans_SILC_Reable=0; %% 1 means the result can be input directly
%         iflag_SLICSeg=0;  %% 1 means the result can be input directly
%         iflag_KmeansSeg=0;  %% 1 means the result can be input directly
        if iflag_Kmeans_SILC_Reable==0 %% =============================search similiar pixels using segmatation results
            if iflag_SLICSeg==0   %%================================================ SLIC segmatation result
                SLIC_CloudArea=CloudAreaSegmatation_SLIC(uint8(J_resample),It_mask,i_Superpixels,Str_SLIC);%%
            else
                SLIC_CloudArea=imread(Str_SLIC);
            end
            if iflag_KmeansSeg==0  %%================================================segmatation result
                KmeansData=uint8(K_means_function(double(J_resample),i_KmeansClasses));%%
                imwrite(uint8(KmeansData), Str_Kmeans);
            else
                KmeansData=imread(Str_Kmeans);
            end
            KmeansData=medfilt2(KmeansData(:,:,1),[5,5], 'symmetric');
            Kmeans_All1=ImgSegmatationRelable_Kmeans(uint8(KmeansData),i_KmeansClasses);   
        %     imwrite(KmeansData,[str_InputPath,'\','J0_Kmeans_OutputLable.tif']);
        %     figure,imshow(40*uint8(Kmeans_All1(:,:,end))) figure,imshow(uint8(OutputLable)) 
            OutputLable=ImgSegmatationRelable_Kmeans_SILC(SLIC_CloudArea,Kmeans_All1,icpend);%icpend 
    %         imwrite(OutputLable,[str_InputPath,'\',str_Jt,'_SILC_Kmeans_OutputLable.png']); 
            WriteMultiBandsImages(uint16(OutputLable),[str_InputPath,'\',str_Jt,'_SILC_Kmeans_OutputLable.png'],16);
        elseif iflag_Kmeans_SILC_Reable==1 %% ===================================
            OutputLable=imread([str_InputPath,'\',str_Jt,'_SILC_Kmeans_OutputLable.png']);
            KmeansData=imread(Str_Kmeans);
            KmeansData=medfilt2(KmeansData(:,:,1),[3,3], 'symmetric');
            Kmeans_All1=ImgSegmatationRelable_Kmeans(KmeansData,i_KmeansClasses); 
        end
        [N_I,N_mask]=GetContralPoints_ChangeDetection_Segmentation_Distance(I,It_mask,J,Jt_mask,OutputLable,Kmeans_All1,str_InputPath,icpend,iwindow_s,[p_StablePixelsSolving{3,2} p_StablePixelsSolving{4,2}]);%icpend
        imwrite(uint8(N_I),[str_InputPath,'\',str_It,'-N_I',str_tif]);
        imwrite(uint8(N_mask),[str_InputPath,'\',str_It,'-N_mask',str_tif]);
        clear SLIC_CloudArea Kmeans_All
        elseif iflag_StablePixels==1
           N_I=imread([str_InputPath,'\',str_It,'-N_I',str_tif]);
           N_mask=imread([str_InputPath,'\',str_It,'-N_mask',str_tif]);
        end
    else
        N_I=I;
        N_mask=It_mask;
    end
    if TarImg==1
        Updated_I=N_I;
        Updated_Imask=N_mask;
    elseif TarImg==2
        Updated_J=N_I;
        Updated_Jmask=N_mask;
    end
end

end
