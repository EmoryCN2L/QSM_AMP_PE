function niftiSetQoffset(filename,offset)
%This function set the sform transformation matrix of a niftifile
%   function niftiSetQoffset(filename,mat)
%       filename,   filename of a nifti image
%       vect,        vector of length 3 
%
%Deqing Qiu. qiudeqiang@gmail.com
%Stanfor University.
%

hdr = niftiLoadHeader(filename);
hdr.hist.qoffset_x = offset(1);
hdr.hist.qoffset_y = offset(2);
hdr.hist.qoffset_z = offset(3);

niftiSaveHeader(hdr,filename);



return