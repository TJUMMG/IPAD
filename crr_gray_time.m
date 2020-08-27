%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2018 by                                          %
%   Jing Liu                                                       %
%   Tianjin University, Tianjin, China                             %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   C. H. Cheng, O. C. Au, C. H. Liu, and K. Y. Yip, "Bit-depth    %
%   expansion by contour region reconstruction," in Proc. IEEE     %
%   Int. Symp. Circuits and Systems (ISCAS), 2009, pp. 944-947.    %
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

function I_crr = crr_gray_time(I_in,H,L)
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

stepratio = downmap./(downmap+upmap);
I_crr = I_in+floor(reshape(stepratio,M,N)*(Te-1));

return;