%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2018 by                                          %
%   Jing Liu                                                       %
%   Tianjin University, Tianjin, China                             %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   P. Wan, G. Cheung, D. Florencio, C. Zhang, and O. C. Au,       %
%   "Image bit-depth enhancement via maximum a posteriori          % 
%   estimation of AC signal," IEEE Trans. Image Processing, vol.   %
%   25, no. 6, pp. 2896-2909, June 2016.                           %
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


function imgACDC = acdc_bits_lambda(imgIn,high,low,blockSize,lambda)
bwidth = 2;
% pre-processing
[ROW,COL,CHN] = size(imgIn);
imgInFloat = imgIn/(2^high-1);
Q = 1/(2^low);
imgInQuant = (floor(imgInFloat/Q)+0.5)*Q;
latest = imgInQuant;

for i = 1:blockSize/2:ROW
    for j = 1:blockSize/2:COL
        endi = min(i+blockSize-1,ROW);
        endj = min(j+blockSize-1,COL);
        block = imgInQuant(i:endi,j:endj,:);
        nDim = (endi-i+1)*(endj-j+1);
        W = generateGraph(block,Q);
        
        D = diag(sum(W,2));
        L = D-W;
        mask = ones(endi-i+1,endj-j+1);
        mask(i+bwidth:endi-i+1-bwidth,j+bwidth:endj-j+1-bwidth)=0;
        bound = find(mask);
        S = zeros(size(bound,1),nDim);
        S(sub2ind([size(bound,1),nDim],(1:size(bound,1))',bound))=1;
        
        % optimization(1 chn)
        for chn = 1:CHN
            blockSingleChn = block(:,:,chn);
            yd = mean(blockSingleChn(:));
            ya = blockSingleChn-yd;
            
            M = [L*Q^2/2/lambda,zeros(nDim,2);...
                zeros(1,nDim),1,-1;zeros(1,nDim),-1,1];
            A = [-eye(nDim),ones(nDim,1),zeros(nDim,1);...
                eye(nDim),zeros(nDim,1),-ones(nDim,1);...
                zeros(1,nDim),-1,1];
            b = [-ya(:);ya(:);Q];
            C = [ones(1,nDim),0,0;S,-0.5*ones(size(bound,1),2)];
            d = [0;S*reshape(latest(i:endi,j:endj,chn),nDim,1)-yd];
            options = optimoptions('quadprog','TolFun',1e-12,'TolX',1e-12);
            [z,fval,exitflag] = quadprog(sparse(M),zeros(nDim+2,1),A,b,C,d,[],[],[],options);
            
            xa = z(1:end-2);
            xd = yd-0.5*(min(xa-ya(:))+max(xa-ya(:)));        
            latest(i:endi,j:endj,chn)=reshape(xa+xd,[endi-i+1,endj-j+1]);
        end
    end
end
imgACDC = latest;

return;
        



