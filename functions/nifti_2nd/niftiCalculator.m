function niftiCalculator(iname,oname,expression,otype)
%This function is a simple nifti image calculator
%function niftiCalculator(iname,oname,expression)
%
%Input:
%   iname,      input name
%   oname,      output name
%   expression, formular expression, e.g. "i1*10" to scale the image by ten
%
%Deqiang Qiu, qiudeqiang@gmail.com
%Emory University, 2013

%%

iImg = niftiLoadNii(iname);
i1 = double(iImg.img)*iImg.hdr.dime.scl_slope + iImg.hdr.dime.scl_inter;
iImg.hdr.dime.scl_slope=1;
iImg.hdr.dime.scl_inter=0;

if(exist('otype','var'))
	iImg.hdr.dime.datatype=otype;
%iImg.hdr
end

eval(sprintf('oImg = %s;',expression));

niftiSaveNii(iImg.hdr,oImg,oname);

return;
