clc, clear
addpath('.\Segmentation_Classification');
str_tif='.tif';
I_str='I0'; 
J_str='J0';
Imask_str='I0_mask';
Jmask_str='J0_mask';
str_cell(1,1)={I_str};
str_cell(1,2)={J_str};
%% =============================== data 1: 1_TrueColorComposition (Reference image J0 is Cloud-free)=============================== 
%%% the Whole Image Cloud Removal by IC_CSC
str_InputPath='..\data\1_TrueColorComposition';
str_cell(1,3)={str_InputPath};
I=imread([str_InputPath,'\',I_str,'_Clouds',str_tif]);
J=imread([str_InputPath,'\',J_str,str_tif]);
I_mask=imread([str_InputPath,'\',Imask_str,str_tif]);
J_mask=I_mask*0; % Reference image J0 is Cloud-free
i_flag=0;
if i_flag==0 %% set "i_flag==0" means without the step of "Stable Pixels Layout", which takes shorter time
    Updated_I=I;
    Updated_Imask=I_mask;
    Updated_J=J;
    Updated_Jmask=J_mask;
elseif i_flag==1 %% set "i_flag==1" means conducting the method with the step of "Stable Pixels Layout", which takes longer time
    p_StablePixelsSolving(1,:)={'the stable pixels are unknown=0 or known=1:',0}; % you can change p_StablePixelsSolving{1,2}=0 if you want to re-solve the stable pixels
    p_StablePixelsSolving(2,:)={'the window size for comparing local similarities:',3}; % you can it to define the window size: eg.3¡Á3
    p_StablePixelsSolving(3,:)={'determine the similarity between I_r(p) and I_r(q):',1.8}; % the threshold to determine the similarity between I_r(p) and I_r(q):(1<t<3)
    p_StablePixelsSolving(4,:)={'determine the consistency between I_r(p) and I_t(p):',0.9}; % the threshold to determine the consistency between I_r(p) and I_t(p):(0.8<t<1)
    p_SLICSeg(1,:)={'the SLIC result is unknown=0 or known=1:',0}; % you can change p_SLICSeg{1,2}=0 if you want to re-do SLIC segmatation
    p_SLICSeg(2,:)={'the size of Superpixels:',15}; % you can change p_SLICSeg{2,2} to define the size of Superpixels in SLIC segmatation
    p_KmeansSeg(1,:)={'the Kmeans result is unknown=0 or known=1:',0}; % you can change p_KmeansSeg{1,2}=0 if you want to re-do Kmeans segmatation
    p_KmeansSeg(2,:)={'the classes of Kmeans:',6}; % you can change p_KmeansSeg{2,2} to define the classes of Kmeans segmatation result
    [Updated_I,Updated_Imask,Updated_J,Updated_Jmask]=Stable_pixels_layout_with_ColorStructure_control(I,I_mask,J,J_mask,str_cell,p_StablePixelsSolving,p_SLICSeg,p_KmeansSeg);
end
[OutputData_II,~]=Whole_Image_Reconstruction_with_Isophote_Constraint(Updated_I,Updated_Imask,Updated_J,Updated_Jmask);
%%% output the result to specified path
imwrite((OutputData_II),[str_InputPath,'\',I_str,'-Result_of_Isophote_',num2str(i_flag),str_tif]);
%% =============================== data 2: 2_FalseColorComposition (both I0 and J0 are Cloudy)=============================== 
%%% the Whole Image Cloud Removal by IC_CSC
str_InputPath='..\data\2_FalseColorComposition';
str_cell(1,3)={str_InputPath};
I=imread([str_InputPath,'\',I_str,'_Clouds',str_tif]);
J=imread([str_InputPath,'\',J_str,'_Clouds',str_tif]);
I_mask=imread([str_InputPath,'\',Imask_str,str_tif]);
J_mask=imread([str_InputPath,'\',Jmask_str,str_tif]);
i_flag=0;
if i_flag==0 %% set "i_flag==0" means without the step of "Stable Pixels Layout", which takes shorter time
    Updated_I=I;
    Updated_Imask=I_mask;
    Updated_J=J;
    Updated_Jmask=J_mask;
elseif i_flag==1 %% set "i_flag==1" means conducting the method with the step of "Stable Pixels Layout", which takes longer time
    p_StablePixelsSolving(1,:)={'the stable pixels are unknown=0 or known=1:',1}; % you can change p_StablePixelsSolving{1,2}=0 if you want to re-solve the stable pixels
    p_StablePixelsSolving(2,:)={'the window size for comparing local similarities:',3}; % you can it to define the window size 
    p_StablePixelsSolving(3,:)={'determine the similarity between I_r(p) and I_r(q):',1.8}; % the threshold to determine the similarity between I_r(p) and I_r(q):(1<t<3)
    p_StablePixelsSolving(4,:)={'determine the consistency between I_r(p) and I_t(p):',0.9}; % the threshold to determine the consistency between I_r(p) and I_t(p):(0.8<t<1)
    p_SLICSeg(1,:)={'the SLIC result is unknown=0 or known=1:',0}; % you can change p_SLICSeg{1,2}=0 if you want to re-do SLIC segmatation
    p_SLICSeg(2,:)={'the size of Superpixels:',15}; % you can change p_SLICSeg{2,2} to define the size of Superpixels in SLIC segmatation
    p_KmeansSeg(1,:)={'the Kmeans result is unknown=0 or known=1:',0}; % you can change p_KmeansSeg{1,2}=0 if you want to re-do Kmeans segmatation
    p_KmeansSeg(2,:)={'the classes of Kmeans:',6}; % you can change p_KmeansSeg{2,2} to define the classes of Kmeans segmatation result
    [Updated_I,Updated_Imask,Updated_J,Updated_Jmask]=Stable_pixels_layout_with_ColorStructure_control(I,I_mask,J,J_mask,str_cell,p_StablePixelsSolving,p_SLICSeg,p_KmeansSeg);
end
[OutputData_II,OutputData_JJ]=Whole_Image_Reconstruction_with_Isophote_Constraint(Updated_I,Updated_Imask,Updated_J,Updated_Jmask);
%%% output the result to specified path
imwrite((OutputData_II),[str_InputPath,'\',I_str,'-Result_of_Isophote_',num2str(i_flag),str_tif]);
imwrite((OutputData_JJ),[str_InputPath,'\',J_str,'-Result_of_Isophote_',num2str(i_flag),str_tif]);
