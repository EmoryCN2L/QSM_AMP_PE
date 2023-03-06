
function niftimathScalarMultiply(fname1,factor,outName);
%function niftimathssub(fname1,fname2,outname);
%this function subtract image2 from image1

[hdr1, filetype1, fileprefix1, machine1] = load_nii_hdr_ut(fname1);


hdrOut = hdr1;
hdrOut.dime.datatype = 16;%float
hdrOut.dime.bitpix = 32;
hdrOut.dime.scl_slop = 1;
hdrOut.dime.scl_interp = 0;
hdrOut.dime.vox_offset = 352;%.nii
hdrOut.hist.magic = 'n+1';

%save header
machineOut = machine1;
outName = imgName(outName);
fidOut = fopen(outName,'wb',machineOut);
assert(fidOut~=-1,sprintf('Fail to open %s for writing',outName));

save_nii_hdr(hdrOut,fidOut);
padHeader(hdrOut,fidOut);

for i = 1:hdrOut.dime.dim(5);
  img1 = load_nii_img(hdr1,filetype1,fileprefix1,machine1,i)*hdr1.dime.scl_slope+hdr1.dime.scl_inter;
  imgOut = img1*factor;
  write_nii_segment(fidOut,imgOut,hdrOut);
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

