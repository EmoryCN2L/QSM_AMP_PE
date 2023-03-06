function o = mypadarray(i, leftsize,rightsize,padval)
%This function pads the array with left size and right size
%function o = mypadarray(i, leftsize,rightsize)
%
%Deqiang Qiu, qiudeqiang@gmail.com
%Stanford University, 2012, Jul

if(~exist('padval','var'))
    padval = 0;
end

o = padarray(padarray(i,leftsize,padval,'pre'),rightsize,padval,'post');

return

