%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2018 by                                          %
%   Jing Liu                                                       %
%   Tianjin University, Tianjin, China                             %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   P. Wan, O. C. Au, K. Tang, Y. Guo, and L. Fang, "From 2D       %
%   extrapolation to 1D interpolation: Content adaptive image      %
%   bit-depth expansion," in Proc. IEEE Int. Conf. Multimedia and  %
%   Expo (ICME), 2012, pp. 170-175.                                %
%                                                                  %
%   This program is free of charge for personal and scientific     %
%   use (with proper citation). The author does NOT give up his    %
%   copyright. Any commercial use is prohibited.                   %
%   YOU ARE USING THIS PROGRAM AT YOUR OWN RISK! THE AUTHOR        %
%   IS NOT RESPONSIBLE FOR ANY DAMAGE OR DATA-LOSS CAUSED BY THE   %
%   USE OF THIS PROGRAM.                                           %
%                                                                  %
%   If you have any questions please contact:                      %
%                                                                  %
%   Jing Liu                                                       %
%   School of Electrical and Information Engineering               %
%   Tianjin University                                             %
%   Tianjin, China                                                 %
%                                                                  %
%                                                                  %
%   email: jliu_tju@tju.edu.cn                                     %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_ca = ca_gray_time(I_in,H,L)
Te = 2^(H-L);
[M,N] = size(I_in);

% 4-neighbor
Iidx = 1:M*N;
UPidx = Iidx-1;
UPidx(mod(UPidx,M)==0) = UPidx(mod(UPidx,M)==0)+1;
DOWNidx = Iidx+1;
DOWNidx(mod(DOWNidx,M)==1) = DOWNidx(mod(DOWNidx,M)==1)-1;
LEFTidx = Iidx-M;
LEFTidx(LEFTidx<1) = LEFTidx(LEFTidx<1)+M;
RIGHTidx = Iidx+M;
RIGHTidx(RIGHTidx>M*N) = RIGHTidx(RIGHTidx>M*N)-M;
% 8-neighbor
ULidx = Iidx-1-M;
ULidx(mod(ULidx,M)==0) = ULidx(mod(ULidx,M)==0)+1+M;
ULidx(ULidx<1) = ULidx(ULidx<1)+1+M;
DLidx = Iidx+1-M;
DLidx(mod(DLidx,M)==1) = DLidx(mod(DLidx,M)==1)-1+M;
DLidx(DLidx<1) = DLidx(DLidx<1)-1+M;
URidx = Iidx+M-1;
URidx(mod(URidx,M)==0) = URidx(mod(URidx,M)==0)-M+1;
URidx(URidx>M*N) = URidx(URidx>M*N)-M+1;
DRidx = Iidx+M+1;
DRidx(mod(DRidx,M)==1) = DRidx(mod(DRidx,M)==1)-M-1;
DRidx(DRidx>M*N) = DRidx(DRidx>M*N)-M-1;

Sup4 = (((I_in(UPidx)-I_in(Iidx))<=Te)&((I_in(UPidx)-I_in(Iidx))>0))|...
    (((I_in(DOWNidx)-I_in(Iidx))<=Te)&((I_in(DOWNidx)-I_in(Iidx))>0))|...
    (((I_in(LEFTidx)-I_in(Iidx))<=Te)&((I_in(LEFTidx)-I_in(Iidx))>0))|...
    (((I_in(RIGHTidx)-I_in(Iidx))<=Te)&((I_in(RIGHTidx)-I_in(Iidx))>0));
Sup8 = (((I_in(ULidx)-I_in(Iidx))<=Te)&((I_in(ULidx)-I_in(Iidx))>0))|...
    (((I_in(DLidx)-I_in(Iidx))<=Te)&((I_in(DLidx)-I_in(Iidx))>0))|...
    (((I_in(URidx)-I_in(Iidx))<=Te)&((I_in(URidx)-I_in(Iidx))>0))|...
    (((I_in(DRidx)-I_in(Iidx))<=Te)&((I_in(DRidx)-I_in(Iidx))>0));
upmap = min(Sup4,Sup8*sqrt(2));
Sdn4 = (((I_in(UPidx)-I_in(Iidx))>=-Te)&((I_in(UPidx)-I_in(Iidx))<0))|...
    (((I_in(DOWNidx)-I_in(Iidx))>=-Te)&((I_in(DOWNidx)-I_in(Iidx))<0))|...
    (((I_in(LEFTidx)-I_in(Iidx))>=-Te)&((I_in(LEFTidx)-I_in(Iidx))<0))|...
    (((I_in(RIGHTidx)-I_in(Iidx))>=-Te)&((I_in(RIGHTidx)-I_in(Iidx))<0));
Sdn8 = (((I_in(ULidx)-I_in(Iidx))>=-Te)&((I_in(ULidx)-I_in(Iidx))<0))|...
    (((I_in(DLidx)-I_in(Iidx))>=-Te)&((I_in(DLidx)-I_in(Iidx))<0))|...
    (((I_in(URidx)-I_in(Iidx))>=-Te)&((I_in(URidx)-I_in(Iidx))<0))|...
    (((I_in(DRidx)-I_in(Iidx))>=-Te)&((I_in(DRidx)-I_in(Iidx))<0));
downmap = min(Sdn4,Sdn8*sqrt(2));


upmap = double(upmap);
upmap(upmap==0) = 100000;
downmap = double(downmap);
downmap(downmap==0) = 100000;

upmap_old = upmap;
while(1)
    tmp = find(I_in(Iidx)==I_in(ULidx));
    upmap(tmp) = min(upmap(tmp),upmap(ULidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(DLidx));
    upmap(tmp) = min(upmap(tmp),upmap(DLidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(URidx));
    upmap(tmp) = min(upmap(tmp),upmap(URidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(DRidx));
    upmap(tmp) = min(upmap(tmp),upmap(DRidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(UPidx));
    upmap(tmp) = min(upmap(tmp),upmap(UPidx(tmp))+1);
    tmp = find(I_in(Iidx)==I_in(DOWNidx));
    upmap(tmp) = min(upmap(tmp),upmap(DOWNidx(tmp))+1);
    tmp = find(I_in(Iidx)==I_in(LEFTidx));
    upmap(tmp) = min(upmap(tmp),upmap(LEFTidx(tmp))+1);
    tmp = find(I_in(Iidx)==I_in(RIGHTidx));
    upmap(tmp) = min(upmap(tmp),upmap(RIGHTidx(tmp))+1);
    
    if all(upmap(:)==upmap_old(:))
        break
    end
    upmap_old = upmap;
end
downmap_old = downmap;
while(1)
    tmp = find(I_in(Iidx)==I_in(ULidx));
    downmap(tmp) = min(downmap(tmp),downmap(ULidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(DLidx));
    downmap(tmp) = min(downmap(tmp),downmap(DLidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(URidx));
    downmap(tmp) = min(downmap(tmp),downmap(URidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(DRidx));
    downmap(tmp) = min(downmap(tmp),downmap(DRidx(tmp))+sqrt(2));
    tmp = find(I_in(Iidx)==I_in(UPidx));
    downmap(tmp) = min(downmap(tmp),downmap(UPidx(tmp))+1);
    tmp = find(I_in(Iidx)==I_in(DOWNidx));
    downmap(tmp) = min(downmap(tmp),downmap(DOWNidx(tmp))+1);
    tmp = find(I_in(Iidx)==I_in(LEFTidx));
    downmap(tmp) = min(downmap(tmp),downmap(LEFTidx(tmp))+1);
    tmp = find(I_in(Iidx)==I_in(RIGHTidx));
    downmap(tmp) = min(downmap(tmp),downmap(RIGHTidx(tmp))+1);
    
    if all(downmap(:)==downmap_old(:))
        break
    end
    downmap_old = downmap;
end
upmap_reshape = reshape(upmap,M,N);
downmap_reshape = reshape(downmap,M,N);

stepratio = downmap./(downmap+upmap);
Thr_upper = 0.99;
Thr_lower = 0.01;
LMin = stepratio>Thr_upper;
LMax = stepratio<Thr_lower;

% denoise1: open to remove tiny LMM regions
se = strel('disk',1);
LMin_opened = imopen(reshape(LMin,M,N),se);
LMax_opened = imopen(reshape(LMax,M,N),se);
% denoise2: removal irregular LMM regions
LMM = LMax_opened+LMin_opened*(-1);
LMM_old = LMM;
L = 1;
while(1)
    for idxi = 1:size(LMM,1)
        for idxj = 1:size(LMM,2)
            if LMM_old(idxi,idxj)~=0
                local = LMM_old(max(idxi-L,1):min(idxi+L,M),...
                    max(idxj-L,1):min(idxj+L,N));
                LMM(idxi,idxj) = 2*double(sum(local(:))>0)-1;
            end
        end
    end
    if all(LMM(:)==LMM_old(:))
        break;
    end
    LMM_old = LMM;
end
% denoise3: close to fill small holes
LMin = LMM<0;
LMax = LMM>0;
LMin_closed = imclose(reshape(LMin,M,N),se);
LMax_closed = imclose(reshape(LMax,M,N),se);
S_Min = bwmorph(LMin_closed,'skel',Inf);
S_Max = bwmorph(LMax_closed,'skel',inf);
S = S_Min|S_Max;

% Virtual Skeleton
% for local minimum
L = double((upmap(Iidx)-upmap(UPidx)>0)&(upmap(Iidx)-upmap(DOWNidx)>0))+...
    double((upmap(Iidx)-upmap(LEFTidx)>0)&(upmap(Iidx)-upmap(RIGHTidx)>0))+...
    double((upmap(Iidx)-upmap(ULidx)>0)&(upmap(Iidx)-upmap(DRidx)>0))+...
    double((upmap(Iidx)-upmap(URidx)>0)&(upmap(Iidx)-upmap(DLidx)>0));
Bd_Min = reshape((L>=2),M,N)&LMM;
L = double((downmap(Iidx)-downmap(UPidx)>0)&(downmap(Iidx)-downmap(DOWNidx)>0))+...
    double((downmap(Iidx)-downmap(LEFTidx)>0)&(downmap(Iidx)-downmap(RIGHTidx)>0))+...
    double((downmap(Iidx)-downmap(ULidx)>0)&(downmap(Iidx)-downmap(DRidx)>0))+...
    double((downmap(Iidx)-downmap(URidx)>0)&(downmap(Iidx)-downmap(DLidx)>0));
Bd_Max = reshape((L>=2),M,N)&LMM;
VS=S|Bd_Min|Bd_Max;

% re-calculate downmap and upmap
upmap_recal = upmap_reshape;
upmap_recal(VS&LMax_closed) = 0;
upmap_recal_old = upmap_recal;
while(1)
    tmp = find((I_in(Iidx)==I_in(ULidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(ULidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(DLidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(DLidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(URidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(URidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(DRidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(DRidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(UPidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(UPidx(tmp))+1);
    tmp = find((I_in(Iidx)==I_in(DOWNidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(DOWNidx(tmp))+1);
    tmp = find((I_in(Iidx)==I_in(LEFTidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(LEFTidx(tmp))+1);
    tmp = find((I_in(Iidx)==I_in(RIGHTidx))&(upmap_reshape(Iidx)==100000));
    upmap_recal(tmp) = min(upmap_recal(tmp),upmap_recal(RIGHTidx(tmp))+1);
    
    if all(upmap_recal(:)==upmap_recal_old(:))
        break
    end
    upmap_recal_old = upmap_recal;
end
upmap_recal(VS&LMax_closed) = 0;

downmap_recal = downmap_reshape;
downmap_recal(VS&LMin_closed) = 0;
downmap_recal_old = downmap_recal;
while(1)
    tmp = find((I_in(Iidx)==I_in(ULidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(ULidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(DLidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(DLidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(URidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(URidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(DRidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(DRidx(tmp))+sqrt(2));
    tmp = find((I_in(Iidx)==I_in(UPidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(UPidx(tmp))+1);
    tmp = find((I_in(Iidx)==I_in(DOWNidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(DOWNidx(tmp))+1);
    tmp = find((I_in(Iidx)==I_in(LEFTidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(LEFTidx(tmp))+1);
    tmp = find((I_in(Iidx)==I_in(RIGHTidx))&(downmap_reshape(Iidx)==100000));
    downmap_recal(tmp) = min(downmap_recal(tmp),downmap_recal(RIGHTidx(tmp))+1);
    
    if all(downmap_recal(:)==downmap_recal_old(:))
        break
    end
    downmap_recal_old = downmap_recal;
end
downmap_recal(VS&LMin_closed) = 0;
stepratio = downmap_recal./(downmap_recal+upmap_recal);

% content adaptive LSBs reconstruction
alpha = 1;
gk = zeros(size(I_in));
gk(LMM==0) = stepratio(LMM==0);
gk(I_in==(2^H-1)-(Te-1)) = cos((1-stepratio(I_in==(2^H-1)-(Te-1))).^alpha);
gk(I_in==0) = 0.5;
gk(LMin_closed==1) = 0.5+0.5*stepratio(LMin_closed==1);
gk(LMax_closed==1) = 0.5*stepratio(LMax_closed==1); 
I_ca2 = I_in+floor(gk*(Te-1));

% bilateral filter
I_ca = I_ca2;
w     = 2;       % bilateral filter half-width
sigma = [2 0.01]; % bilateral filter standard deviations
I_ca_tmp = bfilter2(double(I_ca2)/(2^H-1),w,sigma)*(2^H-1);
I_ca(VS)=I_ca_tmp(VS);
I_ca = round(I_ca);

return;
