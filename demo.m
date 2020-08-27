%  ----------------------------------------------------------------------
%  Executable software for image de-quantization
%
%  ---------------------------------------------------------------------- 
%  Copyright (c) 2018 Jing Liu
%  ----------------------------------------------------------------------
%    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
%  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
%  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%  OTHER DEALINGS IN THE SOFTWARE.
%  ---------------------------------------------------------------------- 
%  Contact: Jing Liu at <jliu_tju@tju.edu.cn> or <jingliu.phd@gmail.com>
%  ----------------------------------------------------------------------
%  Please cite the paper (Jing Liu, Guangtao Zhai, Xiaokang Yang, and 
%  Chang Wen Chen, "IPAD: Intensity Potential for Adaptive De-quantization,
%  " IEEE Trans. Image Process. (TIP), vol.XX, no.X, pp.XXX-XXX, 2018)
%  ----------------------------------------------------------------------

clear; close all; clc;
high = 8; % high bit-depth, change to 16 for 16-bit original image
low = 4; % low bit-depth, smaller than high
maxmemo = 1e8; % limit the memory used
qStep = 2^(high-low); % quantization step
blockSize = 64; % block size for ACDC algorithm
imgName = 'kodim04.png'; % filename for original image
lambda = 1.5; % parameters for ACDC algorithm

ori = imread(imgName);
[M,N,C] = size(ori);
img = double(ori);
img_Q = floor(img/qStep);
% zp
imgIn = img_Q*qStep;
imwrite(uint8(imgIn),sprintf('%s_zp_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(imgIn),sprintf('%s_zp_%d_%d.png',imgName(1:end-4),high,low));
    

% mig
I_mig = img_Q/(2^low-1)*(2^high-1);
imwrite(uint8(I_mig),sprintf('%s_mig_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(I_mig),sprintf('%s_mig_%d_%d.png',imgName(1:end-4),high,low));
  
% br
msbs = uint8(ori)/(2^low);
I_br = imgIn+double(msbs);
imwrite(uint8(I_br),sprintf('%s_br_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(I_br),sprintf('%s_br_%d_%d.png',imgName(1:end-4),high,low));

I_crr = zeros(size(imgIn));
I_ca = zeros(size(imgIn));
I_mrc = zeros(size(imgIn));
I_ipad = zeros(size(imgIn));

% crr
for chn = 1:size(ori,3)
    I_crr(:,:,chn) = crr_gray_time(imgIn(:,:,chn),high,low);
end
imwrite(uint8(I_crr),sprintf('%s_crr_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(I_crr),sprintf('%s_crr_%d_%d.png',imgName(1:end-4),high,low));

% ca
for chn = 1:size(ori,3)
    I_ca(:,:,chn) = ca_gray_time(imgIn(:,:,chn),high,low);
end
imwrite(uint8(I_ca),sprintf('%s_ca_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(I_ca),sprintf('%s_ca_%d_%d.png',imgName(1:end-4),high,low));

% mrc
for chn = 1:size(ori,3)
    I_mrc(:,:,chn) = mrc_gray_bits_fast(imgIn(:,:,chn),high,low);
end
imwrite(uint8(I_mrc),sprintf('%s_mrc_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(I_mrc),sprintf('%s_mrc_%d_%d.png',imgName(1:end-4),high,low));

% acdc
imgACDC = acdc_bits_lambda(imgIn,high,low,blockSize,lambda);
imwrite(uint8(imgACDC*(2^high-1)),sprintf('%s_acdc_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(imgACDC*(2^high-1)),sprintf('%s_acdc_%d_%d.png',imgName(1:end-4),high,low));

% ipad
for chn = 1:size(ori,3)
    I_ipad(:,:,chn) = ipad_gray_time(imgIn(:,:,chn),high,low,maxmemo);
end
imwrite(uint8(I_ipad),sprintf('%s_ipad_%d_%d.png',imgName(1:end-4),high,low));
% imwrite(uint16(I_ipad),sprintf('%s_ipad_%d_%d.png',imgName(1:end-4),high,low));
    