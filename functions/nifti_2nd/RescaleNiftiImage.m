function RescaleNiftiImage(fname,slopeInterp);
%function RescaleNiftiImage(fname,slopeInterp);
%adjust the scale and intercept of an image.



[hdr, filetype, fileprefix, machine] = load_nii_hdr(fname);


switch filetype
    case 1
        ext = '.hdr';
    case 2
        ext = '.nii';
    otherwise
        error('Analyze format not supported');
end;

hdr.dime.scl_slope = slopeInterp(1);
hdr.dime.scl_inter = slopeInterp(2);

fid = fopen([fileprefix,ext],'r+b',machine);
if(fid<0)
    error('Fail to write back file');
end;

fseek(fid,0,'bof');%seek to begin of file

save_nii_hdr(hdr,fid);
fclose(fid);


