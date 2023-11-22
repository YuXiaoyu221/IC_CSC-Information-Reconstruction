function [OutputData_II,OutputData_JJ]=Whole_Image_Reconstruction_with_Isophote_Constraint(Img_I,I_mask,Img_J,J_mask)
%%% I_mask is the cloud detection result of Img_I, where the value 255 indicates the missing area, and the value 0 is the clear area.
%%% J_mask is the cloud detection result of Img_J, where the value 255 indicates the missing area, and the value 0 is the clear area.
for Im=1:2
if (Im==1)
    Img_Target=Img_I;
    Mask_Target=I_mask;
    T0_mask=I_mask; 
    Img_Reference=Img_J; 
end
if (Im==2)
    Img_Target=Img_J;  
    Mask_Target=J_mask;
    T0_mask=J_mask;
    Img_Reference=Img_I; 
end
It_mask=uint8(Mask_Target);
[nh,nw]=size(It_mask);
L_Mask=bwlabel(T0_mask,4);
maxm=max(max(L_Mask));%maxm=1;
OutputData_I=Img_Target;
for nci=1:maxm%maxm
    [ind_h,ind_w]=find(L_Mask==nci);
    J_Mask1=zeros(size(It_mask));
    J_Mask1(L_Mask==nci)=It_mask(L_Mask==nci);
    ssdd=40;
    sy=max(1,min(ind_w)-ssdd);sx=max(1,min(ind_h)-ssdd);ey=min(nw,max(ind_w)+ssdd);ex=min(nh,max(ind_h)+ssdd);
    A_mask=J_Mask1(sx:ex,sy:ey);   
    A0=Img_Target(sx:ex,sy:ey,:);
    B0=Img_Reference(sx:ex,sy:ey,:);
    outi=uint8(AdaptWeightedGradient_E_consistency(A0,A_mask,B0));
    A_Mask10=repmat(A_mask,1,1,3);
    A_Mask100=repmat(J_Mask1,1,1,3);
    A_Mask100(sx:ex,sy:ey,:)=A_Mask10;
    OutputData_I(A_Mask100>0)=outi(A_Mask10>0);
end
if (Im==1)
    OutputData_II=uint8(OutputData_I);
end
if (Im==2)
    OutputData_JJ=uint8(OutputData_I);
end
end

end
