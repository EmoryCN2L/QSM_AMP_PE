function o = myunpadarray(imatrix, leftsize,rightsize)
%This function pads the array with left size and right size
%function o = mypadarray(i, leftsize,rightsize)
%
%Deqiang Qiu, qiudeqiang@gmail.com
%Stanford University, 2012, Jul

if(~exist('padval','var'))
    padval = 0;
end

sz = size(imatrix);
if(numel(sz)<3)
    sz = padarray(sz,3-numel(sz),0,'post');
end

if(numel(sz)>4)
    error('myunpadarray only works for matrix less then 4 dim');
end

sz = toColumnVector(sz);
leftsize = toColumnVector(leftsize);
rightsize = toColumnVector(rightsize);

%pad the size such that they matches the dimension of i
if(numel(leftsize)<numel(sz))
    leftsize = padarray(leftsize,numel(sz)-numel(leftsize),0,'post');
end

if(numel(rightsize)<numel(sz))
    rightsize = padarray(rightsize,numel(sz)-numel(rightsize),0,'post');
end



st = leftsize+1;
ed = sz-rightsize;

if(any(st>ed))
    error('unpad size is larger than matrix size');
end

o = imatrix(st(1):ed(1),st(2):ed(2), st(3):ed(3),:);

return

