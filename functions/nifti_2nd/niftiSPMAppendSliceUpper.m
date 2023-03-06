function niftiSPMAppendSliceUpper(iname,oname,nSlice)
%niftiSPMAppendSlice append slice at the end of a spm volume without
%changing it's location
%function niftiSPMAppendSlice(iname,oname,nSlice)
%iname, input name
%oname, output name
%nSlice, number of slices to append

hdr = spm_vol(iname);

image = zeros(hdr.dim(1),hdr.dim(2),hdr.dim(3)+nSlice);

image(:,:,1:hdr.dim(3)) = spm_read_vols(hdr);


hdr.dim(3) = hdr.dim(3)+nSlice;

hdr.fname = oname;

spm_write_vol(hdr,image);
return;