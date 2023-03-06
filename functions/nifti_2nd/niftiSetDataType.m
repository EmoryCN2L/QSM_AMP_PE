function ohdr = niftiSetDataType(ihdr,t)
%set nifti data type of the header hdr
%hdr = niftiSetDataType(hdr,t)
%ohdr, output header with datatype properly set
%ihdr, input header
%t, data type: can be one of
% 'uint8','int16','int32','float32','float64','int8','uint16','uint32'

ohdr = ihdr;
ohdr.dime.datatype = niftiType(t);
ohdr.dime.bitpix = niftiType(t,'bits');

return;
