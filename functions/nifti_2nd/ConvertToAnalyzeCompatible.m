function ConvertToAnalyzeCompatible(iname,oname);
%change the storage order to be analyze compatible
%function ConvertToAnalyzeCompatible(iname,oname);
%

[ihdr,ift,ifp,imachine] = load_nii_hdr(iname);

if(ihdr.dime.pixdim(1)==-1)
  %already in analyze storage order copy file
  copyfile(iname,oname);
  return;
end;

ohdr = ihdr;
ohdr.dime.pixdim(1) = -1;%flip qform
ohdr.hist.srow_x(1) = -1*ohdr.hist.srow_x(1);%flip qform

ofid = fopen(oname,'wb');
save_nii_hdr(ohdr,ofid);

niftiPadHeader(ohdr,ofid);

for i = 1:ohdr.dime.dim(5);
  img = load_nii_img(ihdr,ift,ifp,imachine,i);
  img = flipdim(double(img),1);
  write_nii_segment(ofid, img, ohdr);
  disp([num2str(i),'completed']);
end;

fclose(ofid);

return;
