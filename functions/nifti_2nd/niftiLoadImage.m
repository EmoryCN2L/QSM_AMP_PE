%function [img,hdr] = niftiLoadImage(hdr,filetype,fileprefix,machine,img_idx,dim5_idx,dim6_idx,dim7_idx,old_RGB)
%internal function
%This function 

%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

function [img,hdr] = niftiLoadImage(varargin)
    if(numel(varargin)==1)
        [hdr,filetype, fileprefix,machine] = niftiLoadHeader(varargin{1});
        [img,hdr] = load_nii_img(hdr,filetype,fileprefix,machine);
    else
        [img,hdr] = load_nii_img(varargin{:});
    end
return						% read_image

