function niftiCopyQFormSpace(orgname,filename)
%This function copy the qform transformation matrix of a niftifile to
%another
%   function niftiSetSFormSpace(filename,mat)
%       filename,   filename of a nifti image
%       mat,        4x4 matrix
%
%Deqing Qiu. qiudeqiang@gmail.com
%Stanfor University.
%
[hdrOrg,typeOrg] = niftiLoadHeader(orgname);

[hdr,type] = niftiLoadHeader(filename);
hdr.hist.qform_code = hdrOrg.hist.qform_code;
hdr.hist.quatern_b = hdrOrg.hist.quatern_b;
hdr.hist.quatern_c = hdrOrg.hist.quatern_c;
hdr.hist.quatern_d = hdrOrg.hist.quatern_d;
hdr.hist.qoffset_x = hdrOrg.hist.qoffset_x;
hdr.hist.qoffset_y = hdrOrg.hist.qoffset_y;
hdr.hist.qoffset_z = hdrOrg.hist.qoffset_z;
hdr.dime.pixdim(1:4) = hdrOrg.dime.pixdim(1:4);

if(type==2)
    img = niftiLoadImage(filename);
    niftiSaveNii(hdr,img,filename);
else%nifti pair
    niftiSaveHeader(hdr,filename);
end



return