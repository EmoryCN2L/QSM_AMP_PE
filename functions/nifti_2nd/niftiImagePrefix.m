function [prefix,tp] = niftiPrefixAndTypeFromName(fileName);
%0; analyze
%1: nii/.hdr .img;
%2: nii/.nii

%if file extension is not one of ".hdr", ".img" or '.nii'
%.nii is assumed

[p n e] = fileparts(fileName);

if(strcmpi(e,'.hdr')||strcmpi(e,'.img'))
    prefix = fullfile(p,n);
    tp = 1;
elseif(strcmpi(e,'.nii'))
    prefix = fullfile(p,n);
    tp = 2;
else
    prefix = fileName;
    tp = 2;
end;

return;