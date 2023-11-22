function OutputLable=ImgSegmatationRelable_Kmeans_SILC(SLIC_CloudArea,Kmeans_All,Num_clouds)
%%% based on the result of SILC and Kmeans, relable
%%% SLIC_CloudArea: unit16 type; size(~,3) indicates the number of cloud areas, band i is the result of cloud i
%%% Kmeans_All: unit8 type; size(~,3)=class_num+1, band_i is the result of class i
%%% OutputLable: unit16 type; size(~,3) indicates the number of cloud areas, band i is the result of SILC¡ÉK_means
Kmeans_All=Kmeans_All(:,:,end);
[nm,nn,nb]=size(SLIC_CloudArea);nb=min(Num_clouds,nb); %%nb: num of cloud areas
OutputLable=uint16(zeros(nm,nn,size(SLIC_CloudArea,3)));
OutputLable(:,:,nb)=SLIC_CloudArea(:,:,nb);
for ibb=1:nb
    OutputLable_ibb=OutputLable(:,:,ibb);%%fill it class-by-class
    ClassNumbers=0;%% the num of small regions in current cloud
    SLIC_Cloud_i=SLIC_CloudArea(:,:,ibb);
    icl_n=max(max(SLIC_Cloud_i));%
    for icl_i=1:icl_n
        SLIC_Cloud_I_pi=uint8(zeros(nm,nn));
        %%%the num of Kmeans classes in current window
        SLIC_Cloud_I_pi(SLIC_Cloud_i==icl_i)=Kmeans_All(SLIC_Cloud_i==icl_i);       
        %%%pixelValue£¬0~255£¨x-£©£»count the num£¨y-£©
        [count,~]=imhist(SLIC_Cloud_I_pi(SLIC_Cloud_I_pi>0));
        indexi=find(count>0);%%the num of Kmeans classes: index=[3,6],indicating class 2 and class 5
        for ikk=1:length(indexi) %%
            temp_ik=zeros(nm,nn);
            temp_ik(SLIC_Cloud_I_pi==indexi(ikk)-1)=1;
            La_Mask=bwlabel(temp_ik,4);Num_C=max(max(La_Mask));%%
            for Num_C_ii=1:Num_C
                ClassNumbers=ClassNumbers+1;
                OutputLable_ibb(La_Mask==Num_C_ii)=ClassNumbers;
            end
        end
    end
    OutputLable(:,:,ibb)=OutputLable_ibb; %figure,imshow(uint8(OutputLable_ibb))
end
end