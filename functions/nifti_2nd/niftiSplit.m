function niftisplit(iname, oprefix);
%function niftisplit(iname, oname);
%split nifti image in the 4th dimension

[ihdr,ifiletype,iprefix,imachine] = load_nii_hdr_ut(iname);

switch ifiletype
  case 1
    extImg = '.img';
  case 2;
    extImg = '.nii';
  otherwise 
    error('Convertional Analyze format not supported');
end;

fidImgIn = fopen([iprefix,extImg],'rb',imachine);
blocksize = prod(ihdr.dime.dim(2:4))*ihdr.dime.bitpix/8;

ohdr =ihdr;
ohdr.dime.vox_offset = 352;
ohdr.dime.dim(5) = 1;


status = fseek(fidImgIn,ihdr.dime.vox_offset,'bof');
assert(status==0,'Fail to move file pointer');

for i = 1:ihdr.dime.dim(5);
  %open a .nii file for writing
  fidImgOut = fopen(sprintf('%s%04d.nii',oprefix,i-1),'wb',imachine);
  assert(fidImgOut~=-1, 'Fail to open image for writing');
  
  save_nii_hdr(ohdr,fidImgOut);
  curPos = ftell(fidImgOut);
  
  bytewrote = fwrite(fidImgOut,zeros(ohdr.dime.vox_offset-curPos,1),'uchar');%pad header

  %read data and write out
  data = fread(fidImgIn,blocksize,'uchar=>uchar');
  bytewrote = fwrite(fidImgOut,data,'uchar');
  assert(bytewrote==blocksize,'Incomplete image wrote');
  
  fclose(fidImgOut);  
  disp([num2str(i), 'images completed']);
end;

fclose(fidImgIn);

