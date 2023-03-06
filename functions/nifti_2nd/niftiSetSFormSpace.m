function niftiSetSFormSpace(filename,mat)
%This function set the sform transformation matrix of a niftifile
%   function niftiSetSFormSpace(filename,mat)
%       filename,   filename of a nifti image
%       mat,        4x4 matrix
%
%Deqing Qiu. qiudeqiang@gmail.com
%Stanfor University.
%

[hdr,type] = niftiLoadHeader(filename);
hdr.hist.srow_x = mat(1,:);
hdr.hist.srow_y = mat(2,:);
hdr.hist.srow_z = mat(3,:);
hdr.hist.sform_code = 1;
if(type==2)
    img = niftiLoadImage(filename);
    niftiSaveNii(hdr,img,filename);
else%nifti pair
    niftiSaveHeader(hdr,filename);
end



return