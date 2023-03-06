function hdr = niftiSetFileType(hdr,type)
%Purpose: set the magic field according to image type
%Syntax: function hdr = niftiSetFileType(hdr,type)
%Para:
%hdr, file handle
%Para:
%type, numerical or file name
%0: analyze
%1: nii/.hdr .img
%2: nii/.nii

if(ischar(type))
    [prefix,type] = niftiImagePrefix(type);
end;

switch type;
    case 0:
        hdr.hist.magic = '';
    case 1:
        hdr.hist.magic = 'ni1';
    case 2:
        hdr.hist.magic = 'n+1';
end;
return;
        

