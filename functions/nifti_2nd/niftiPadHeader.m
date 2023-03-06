function niftiPadHeader(hdr,fid);
%pad from current file position to the position of hdr.dime.vox_offset
curPos = ftell(fid);
if(curPos ==hdr.dime.vox_offset)
  return;
elseif(curPos>hdr.dime.vox_offset)
  error('Larger header written than specified');
else
  fwrite(fid,zeros(hdr.dime.vox_offset-curPos,1),'uchar');
end;

return;
