function SLIC_CloudArea=CloudAreaSegmatation_SLIC(ImageData,I_mask,ss,Str_Output)
%%%segmatation for Image I with cloud mask "I_mask"
%%% ss is the size of superpixel
%%% Str_Output is the output path
[nh,nw]=size(I_mask);
L_Mask=bwlabel(I_mask,4);
maxm=max(max(L_Mask));
SLIC_CloudArea=zeros(nh,nw,maxm);
for nci=1:maxm
    [ind_h,ind_w]=find(L_Mask==nci);
    J_Mask1=zeros(size(I_mask));
    J_Mask1(L_Mask==nci)=1;
    ssdd=5;
    sy=max(1,min(ind_w)-ssdd);sx=max(1,min(ind_h)-ssdd);ey=min(nw,max(ind_w)+ssdd);ex=min(nh,max(ind_h)+ssdd);
    A_mask=uint8(J_Mask1(sx:ex,sy:ey));
    B0=ImageData(sx:ex,sy:ey,:);
    
    sh=floor(size(B0.*A_mask,1)/ss);
    sw=floor(size(B0.*A_mask,2)/ss);
    if (sh*sw==0)
        outi0=double(A_mask);
    else
        outi0=SLIC_function(B0.*A_mask,ss);
    end
    outi=outi0;
    %%%re-lable
    outi=outi.*double(A_mask);
    re_lable=0;
    for ppi=1:max(max(outi))
        a=find(outi==ppi);
        if (~isempty(a))
            re_lable=re_lable+1;
            outi(a)=-1*re_lable;
        end
    end
    SLIC_CloudArea(sx:ex,sy:ey,nci)=-1*outi;
end
% imwrite(uint16(SLIC_CloudArea),Str_Output,'bitdepth',16);
WriteMultiBandsImages(uint16(SLIC_CloudArea),Str_Output,16);
end