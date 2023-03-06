function niftiVOI(fname,outName,volumes)
%function niftiVOI(fname1,volumes)
[hdr1, filetype1, fileprefix1, machine1] = load_nii_hdr_ut(fname);

hdrOut = hdr1;
hdrOut.dime.dim(5) = numel(volumes);

fidOut = fopen(imgName(fname),'wb',machine1);
save_nii_hdr(hdrOut,fidOut);
padHeader(hdrOut,fidOut);

for i = 1:numel(volumes);
img = load_nii_img(hdr1,filetype1,fileprefix1,machine1,volumes(i));
write_nii_segment(fidOut,img,hdrOut);
disp(sprintf('%d finished',i));
end;


return;





%=============================================
%append .nii to the end of file name as neccesary
function o = imgName(fname)
  [p n e] = fileparts(fname);
  if(~strcmpi(e,'.nii'))
    o = [fname,'.nii'];
  else
    o = fname;
  end;
return;


%=============================================
function padHeader(hdr,fid)
%
curPos = ftell(fid);
if(curPos<hdr.dime.vox_offset)
  fwrite(fid,zeros(hdr.dime.vox_offset-curPos,1),'uchar');
end;

return;