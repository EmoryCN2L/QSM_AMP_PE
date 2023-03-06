function niftiAppendSlice(iname,oname)
%This function append slice at the bottom of the original image but retain
%the coordinate
%%function niftiApendSlice(iname,ns,oname);
%iname,input name of nifti image
%ns, number of slice
%oname, outputname of nifti image
%
%@Copyright, Deqiang Qiu (qiudeqiang@gmail.com), Stanford University


ihdr = niftiLoadHeader(iname);
