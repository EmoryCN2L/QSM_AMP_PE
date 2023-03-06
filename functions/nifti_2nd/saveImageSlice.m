function saveImageSlice(img,name);


hdr = spm_vol;
hdr(1).dim = [size(img,1),size(img,2),size(img,3)];
hdr.mat = eye(4);

hdr.dt = [spm_type('float32'),0];
hdr.pinfo = [1,0,0]';
hdr.fname = name;

spm_write_vol(hdr,img);

return;