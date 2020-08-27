function W = generateGraph(block, Q)

[blockM,blockN,~]=size(block);
nDim = blockM*blockN;
% sparse version
W = zeros(nDim,nDim);
mask = reshape(1:nDim,[blockM,blockN]);
list = reshape(mask(1:blockM-1,1:blockN-1),1,(blockM-1)*(blockN-1));
for i = 1:size(list,2)
    W(list(i),list(i)+1)=max(abs([block(list(i)),block(list(i)+nDim),block(list(i)+2*nDim)]-...
        [block(list(i)+1),block(list(i)+1+nDim),block(list(i)+1+2*nDim)]))<=Q;
    W(list(i),list(i)+blockM)=max(abs([block(list(i)),block(list(i)+nDim),block(list(i)+2*nDim)]-...
        [block(list(i)+blockM),block(list(i)+blockM+nDim),block(list(i)+blockM+2*nDim)]))<=Q;
end
W = W+W';
W = sparse(W);
return;

