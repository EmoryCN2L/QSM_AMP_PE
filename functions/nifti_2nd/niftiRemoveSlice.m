function niftiRemoveSlice(iname,oname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name: niftiRemoveslice
% 
% Desc: Remove slices
% 
% Inputs:  
%   iname,      input file name
%   oname,      output file name
% 
% Outputs: 
% 
% Created: 30-Oct-2011 13:59:16
% 
% Author:  Deqiang Qiu, qiudeqiang@gmail.com
% Stanford University
% All Rights reserved.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hdr = spm_vol(iname);
img = spm_read_vols(hdr);

ohdr = hdr;
ohdr.fname = oname;

img(:,:,1) = 0;
img(:,:,end) = 0;

spm_write_vol(ohdr,img);