function niftiSaveNii(hdr,img,filename,filetype)
%niftiSaveNii(hdr,img,filename,filetype)
%Save Nifti Image,
%if filetype is specified then the designated file type will be used
%irrespect of the extension of file name;
%Convention: .hdr .img /nifti image pair; .nii single niftiImage;
%
%specification of filetype will overwrite the extension
%
if(exist('filetype','var')==0)
    filetype =[];
end;

[p, n, e] = fileparts(filename);


[prefix,tp] = niftiPrefixAndTypeFromName(filename);

if (isempty(filetype))
    filetype = tp;
end;


niftiSaveHeader(hdr,filename,filetype);
fid = niftiOpenImageFileForWrite(filename,hdr,filetype);
niftiSaveImageSegment(fid,img,hdr);
fclose(fid);
return;
