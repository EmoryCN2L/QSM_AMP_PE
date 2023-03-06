function imgMergeFast(outname, fname1,fname2);
%merge the images together in t direction

[hdr1, filetype1, fileprefix1, machine1] = load_nii_hdr(fname1);
[hdr2, filetype2, fileprefix2, machine2] = load_nii_hdr(fname2);

if(any(hdr1.dime.dim(2:4)~=hdr2.dime.dim(2:4))...
    ||hdr1.dime.datatype~=hdr2.dime.datatype...
    ||hdr1.dime.scl_slope~=hdr2.dime.scl_slope...
    ||hdr1.dime.scl_inter~=hdr2.dime.scl_inter...
    ||hdr1.hist.qform_code~=hdr2.hist.qform_code...
    ||hdr1.hist.sform_code~=hdr2.hist.sform_code...
    ||hdr1.hist.quatern_b~=hdr2.hist.quatern_b...
    ||hdr1.hist.quatern_c~=hdr2.hist.quatern_c...
    ||hdr1.hist.quatern_d~=hdr2.hist.quatern_d...
    ||hdr1.hist.qoffset_x~=hdr2.hist.qoffset_x...
    ||hdr1.hist.qoffset_y~=hdr2.hist.qoffset_y...
    ||hdr1.hist.qoffset_z~=hdr2.hist.qoffset_z...
    ||any(hdr1.hist.srow_x~=hdr2.hist.srow_x)...
    ||any(hdr1.hist.srow_y~=hdr2.hist.srow_y)...
    ||any(hdr1.hist.srow_z~=hdr2.hist.srow_z))
  error('Inconsistent image dimensions or storage type');
end;


%save output header, assume nii format
hdrOut = hdr1;
hdrOut.dime.dim(5) = hdr1.dime.dim(5)+hdr2.dime.dim(5);
ofid = fopen(outname,'wb',machine1);
save_nii_hdr(hdrOut,ofid);
curPos = ftell(ofid);

%pad file
if(curPos<hdrOut.dime.vox_offset)
  fwrite(ofid,zeros((hdrOut.dime.vox_offset-curPos),1),'uchar');
end;


%copy the data from first volume
fid1 = fopen(fname1,'rb',machine1);
fseek(fid1,hdr1.dime.vox_offset,'bof');

blockSize = prod(hdr1.dime.dim(2:4))*hdr1.dime.bitpix/8;

for i = 1:hdr1.dime.dim(5);%loop over volumes
  data = fread(fid1,blockSize,'uchar=>uchar');
  byteWrote = fwrite(ofid,data,'uchar');
  if(byteWrote~=blockSize)
    error('Fail to write file');
  end;
  disp([num2str(i/hdrOut.dime.dim(5)*100),'%completed']);
end;

fclose(fid1);

%copy data from the second volume
fid2 = fopen(fname2,'rb',machine2);
fseek(fid2,hdr2.dime.vox_offset,'bof');

for i = 1:hdr2.dime.dim(5);%loop over volumes
  data = fread(fid2,blockSize,'uchar=>uchar');
  byteWrote = fwrite(ofid,data,'uchar');
  if(byteWrote~=blockSize)
    error('Fail to write file');
  end;
  disp([num2str((i+hdr1.dime.dim(5))/hdrOut.dime.dim(5)*100),'%completed']);
end;

fclose(fid2);

fclose(ofid);



