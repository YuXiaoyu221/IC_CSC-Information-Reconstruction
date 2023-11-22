function OutputData2=AdaptWeightedGradient_E_consistency(inputData00,DataMask00, inputData10)
%%%The reciprocal of the square of the gradient is defined as the weight
%%%============= same gradient with same weight £¡£¡£¡======================
[m,n,bs]=size(inputData00);
inputData0=zeros(m+2,n+2,bs)-1;inputData0(2:end-1,2:end-1,:)=inputData00;
inputData1=zeros(m+2,n+2,bs)-1;inputData1(2:end-1,2:end-1,:)=inputData10;
DataMask0=imfilter(DataMask00, [1,1,1;1,1,1;1,1,1]);%mask expansion
DataMask0(DataMask0>0)=1;
DataMask=zeros(m+2,n+2,1)-1;DataMask(2:end-1,2:end-1,:)=DataMask0;
DataMask0=zeros(m+2,n+2,1)-1;DataMask0(2:end-1,2:end-1,:)=DataMask00;
[Array_H,Array_W]=find(DataMask>0);%

[m,n,bs]=size(inputData0);
% DataMask=imfilter(DataMask0, [1,1,1;1,1,1;1,1,1]);%
% DataMask(DataMask>0)=1;
% [Array_H,Array_W]=find(DataMask>0);%
%%% matrix DiagA
%%% the value of boundary pixels in DiagA is 1
%%%creat sparse matrix
% DiagA=eye(length(Array_H));
DiagA_i=zeros(1,length(Array_H));
DiagA_j=zeros(1,length(Array_H));
DiagA_v=ones(1,length(Array_H));
ArrayWH_index=zeros(max(Array_H),max(Array_W));
for i=1:length(Array_W)
    ArrayWH_index(Array_H(i),Array_W(i))=i;
end
OutputData2=inputData0;
for ib=1:bs
    inputData1_ib=double(inputData1(:,:,ib));
    %the gradient in different directions
    d_inputData1=zeros(8,m*n);%
    t=imfilter(double(inputData1_ib), [1,0,0;0,-1,0;0,0,0]);%% left-up
    d_inputData1(1,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,1,0;0,-1,0;0,0,0]);%% up
    d_inputData1(2,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,0,1;0,-1,0;0,0,0]);%% right-up
    d_inputData1(3,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,0,0;1,-1,0;0,0,0]);%% left
    d_inputData1(4,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,0,0;0,-1,1;0,0,0]);%% right
    d_inputData1(5,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,0,0;0,-1,0;1,0,0]);%% left-bottom
    d_inputData1(6,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,0,0;0,-1,0;0,1,0]);%% bottom
    d_inputData1(7,:)=reshape(t,1,m*n);
    t=imfilter(double(inputData1_ib), [0,0,0;0,-1,0;0,0,1]);%% right-bottom
    d_inputData1(8,:)=reshape(t,1,m*n);
%     d_inputData1(1,:)=0;d_inputData1(3,:)=0;d_inputData1(6,:)=0;d_inputData1(8,:)=0;
    % The reciprocal of the square of the gradient is defined as the weight
      d_inputData1(d_inputData1==0)=0.1;%%
  % d_inputData1=d_inputData1+0.01;
    d_inputData1=abs(1./d_inputData1./d_inputData1);% 
    GradientWeights=zeros(m,n,8);
    for ii=1:8
        GradientWeights(:,:,ii)=reshape(d_inputData1(ii,:),m,n);
    end
%     GradientWeights=ones(size(GradientWeights))*0.25;
    clear d_inputData1 d_inputData1_t sumt %% KK
    OutputData=zeros(m,n);
    num_nozeros=1;%the num of non-zero values
for i=1:length(Array_W)
    if ( DataMask0(Array_H(i),Array_W(i))>0)%% unknown pixels
%         if ( DataMask(Array_H(i)-1,Array_W(i))>0 && DataMask(Array_H(i)+1,Array_W(i))>0 && DataMask(Array_H(i),Array_W(i)-1)>0 && DataMask(Array_H(i),Array_W(i)+1)>0 && ...
%        DataMask(Array_H(i)-1,Array_W(i)-1)>0 && DataMask(Array_H(i)+1,Array_W(i)-1)>0 && DataMask(Array_H(i)-1,Array_W(i)+1)>0 && DataMask(Array_H(i)+1,Array_W(i)+1)>0 ) 
        if ( DataMask(Array_H(i)-1,Array_W(i))>0 && DataMask(Array_H(i)+1,Array_W(i))>0 && DataMask(Array_H(i),Array_W(i)-1)>0 && DataMask(Array_H(i),Array_W(i)+1)>0 && ...
       DataMask(Array_H(i)-1,Array_W(i)-1)>0 && DataMask(Array_H(i)+1,Array_W(i)-1)>0 && DataMask(Array_H(i)-1,Array_W(i)+1)>0 && DataMask(Array_H(i)+1,Array_W(i)+1)>0)            
   %solve A in Ax=b, Diag=A
            x1=GradientWeights(Array_H(i),Array_W(i),1);
            x2=GradientWeights(Array_H(i),Array_W(i),2);
            x3=GradientWeights(Array_H(i),Array_W(i),3);
            x4=GradientWeights(Array_H(i),Array_W(i),4);  
            x5=GradientWeights(Array_H(i),Array_W(i),5);
            x6=GradientWeights(Array_H(i),Array_W(i),6);
            x7=GradientWeights(Array_H(i),Array_W(i),7);
            x8=GradientWeights(Array_H(i),Array_W(i),8); 
            x0=(x1+x2+x3+x4+x5+x6+x7+x8);
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i; DiagA_v(num_nozeros)=-x0; num_nozeros=num_nozeros+1;   
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)-1); DiagA_v(num_nozeros)=x1; num_nozeros=num_nozeros+1;
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i-1; DiagA_v(num_nozeros)=x2; num_nozeros=num_nozeros+1;
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)+1); DiagA_v(num_nozeros)=x3; num_nozeros=num_nozeros+1;
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)-1); DiagA_v(num_nozeros)=x4; num_nozeros=num_nozeros+1;
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)+1); DiagA_v(num_nozeros)=x5; num_nozeros=num_nozeros+1;
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)-1); DiagA_v(num_nozeros)=x6; num_nozeros=num_nozeros+1;
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i+1; DiagA_v(num_nozeros)=x7; num_nozeros=num_nozeros+1;  
            DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)+1); DiagA_v(num_nozeros)=x8; num_nozeros=num_nozeros+1;        
            %solve b in Ax=b, OutputData=b
            t=x1*inputData1_ib(Array_H(i)-1,Array_W(i)-1)+x2*inputData1_ib(Array_H(i)-1,Array_W(i))+x3*inputData1_ib(Array_H(i)-1,Array_W(i)+1)+...
                x4*inputData1_ib(Array_H(i),Array_W(i)-1)+x5*inputData1_ib(Array_H(i),Array_W(i)+1)+...
                x6*inputData1_ib(Array_H(i)+1,Array_W(i)-1)+x7*inputData1_ib(Array_H(i)+1,Array_W(i))+x8*inputData1_ib(Array_H(i)+1,Array_W(i)+1)-...
                x0*inputData1_ib(Array_H(i),Array_W(i));
            OutputData(Array_H(i),Array_W(i))=t;
        else %%±ß½çµã
              s=0;ls=0;
              if(DataMask(Array_H(i)-1,Array_W(i)-1)>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)-1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i)-1,Array_W(i)-1,:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i)-1,Array_W(i))>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i-1; DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i)-1,Array_W(i),:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i)-1,Array_W(i)+1)>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)-1,Array_W(i)+1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i)-1,Array_W(i)+1,:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i),Array_W(i)-1)>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)-1); DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1; 
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i),Array_W(i)-1,:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i),Array_W(i)+1)>0)   
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i),Array_W(i)+1); DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i),Array_W(i)+1,:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i)+1,Array_W(i)-1)>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)-1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i)+1,Array_W(i)-1,:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i)+1,Array_W(i))>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i+1; DiagA_v(num_nozeros)=5/16; num_nozeros=num_nozeros+1;
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i)+1,Array_W(i),:)*DiagA_v(num_nozeros-1);
              end
              if(DataMask(Array_H(i)+1,Array_W(i)+1)>0)
                  DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=ArrayWH_index(Array_H(i)+1,Array_W(i)+1); DiagA_v(num_nozeros)=-1/16; num_nozeros=num_nozeros+1;  
                  s=s+DiagA_v(num_nozeros-1);
                  ls=ls+inputData1_ib(Array_H(i)+1,Array_W(i)+1,:)*DiagA_v(num_nozeros-1);
              end
              DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i; DiagA_v(num_nozeros)=-s; num_nozeros=num_nozeros+1;
              OutputData(Array_H(i),Array_W(i),:)=ls-s*inputData1_ib(Array_H(i),Array_W(i),:);
        end
    else
        DiagA_i(num_nozeros)=i; DiagA_j(num_nozeros)=i; DiagA_v(num_nozeros)=1; num_nozeros=num_nozeros+1;
        OutputData(Array_H(i),Array_W(i))=inputData0(Array_H(i),Array_W(i),ib);
    end
end 
    DiagA=sparse(DiagA_i,DiagA_j,DiagA_v);
    clear inputData1_ib GradientWeights DiagA_i DiagA_j DiagA_v;
    OutputData1=OutputData(DataMask>0);
    OutputData1=double(OutputData1);
%     setup = struct('type','ilutp','droptol',1e-6);
%     [L,U] = ilu((DiagA),setup);i
 %   OutputData1=bicgstab(DiagA,OutputData1,1e-4,100,L,U);
    OutputData1=DiagA\OutputData1;
    OutputData=double(inputData0(:,:,ib));
    OutputData(DataMask>0)=abs(OutputData1);
    OutputData2(:,:,ib)=OutputData;
%     OutputData2=zeros(size(DiagA));OutputData2(DiagA~=0)=DiagA(DiagA~=0);
%     inputData1(663:665,1559:1561,:)  OutputData(663:665,1559:1561,:)
end
OutputData2=OutputData2(2:end-1,2:end-1,:);
OutputData2=uint8(OutputData2);
end
% (a(1,1)-a(2,2))/100+(a(1,2)-a(2,2))/100+(a(1,3)-a(2,2))/1+(a(2,1)-a(2,2))/16+(a(2,3)-a(2,2))/100+(a(3,1)-a(2,2))/49+(a(3,2)-a(2,2))/9+(a(3,3)-a(2,2))/100