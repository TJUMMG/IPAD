%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2018 by                                          %
%   Jing Liu                                                       %
%   Tianjin University, Tianjin, China                             %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   G. Mittal, V. Jakhetiya, S. P. Jaiswal, O. C. Au, A. K. Tiwari,%
%   and D. Wei, "Bit-depth expansion using minimum risk based      %
%   classification," in Proc. IEEE Visual Communications and Image %
%   Processing (VCIP), 2012, pp. 1-5.                              %
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

function I_mrc=mrc_gray_bits_fast(I_in,High,Low)
% only for gray image
I_in = double(I_in);
[X,Y] = size(I_in);
I_zp = I_in;
I_mig = round(floor(I_in/(2^(High-Low)))*(2^High-1)/(2^Low-1));
I_mrc = I_mig;

mean_image = round((I_mig(1:end-2,2:end-1) + I_mig(3:end,2:end-1) ...
    + I_mig(2:end-1,1:end-2) + I_mig(2:end-1,3:end))/4); 
Error_image = I_mig(2:end-1,2:end-1)-mean_image;
error_image = reshape(Error_image,1,(X-2)*(Y-2));
mean_image = reshape(mean_image,1,(X-2)*(Y-2));
zp_image = reshape(I_zp(2:X-1,2:Y-1),1,(X-2)*(Y-2));
XX = -(2^High-1):(2^High-1);
kkk = hist(error_image,XX);

pm = zeros(2^(High-Low),(X-2)*(Y-2));
for k = 1:2^(High-Low)
    pm(k,:) = kkk(2^High+(zp_image+k-1-mean_image))';
end
loss = zeros(2^(High-Low),2^(High-Low));
for k=1:2^(High-Low)
    loss(k,:) = ([1:2^(High-Low)]-k).^2;
end

R = loss*pm;
[~,minIdx] = min(R,[],1);
I_mrc(2:end-1,2:end-1) = I_zp(2:end-1,2:end-1)+reshape(minIdx,[X-2,Y-2])-1;
I_mrc = round(I_mrc);
