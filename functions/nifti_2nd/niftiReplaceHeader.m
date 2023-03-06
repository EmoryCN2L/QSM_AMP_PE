function niftiReplaceHeader(hdr,fileName)
%
[prefix,type] = niftiImagePrefix(fileName);

switch type;
    case 1
        hdrName = [prefix,'.hdr'];
        fid = fopen(hdrName,'wb');
        assert(fid~=-1,sprintf('Fail to open file %s for writing',hdrName));
        save_nii_hdr(hdr,fid);
        fclose(fid);
    case 2
        hdrName = [prefix,'.nii'];
        fid = fopen(hdrName,'r+b');
        assert(fid~=-1,sprintf('Fail to open file %s for writing',hdrName));
        fseek(fid,0,'bof');
        pos = ftell(fid);
        assert(pos==0,'Fail to seek to beginning of file');
        save_nii_hdr(hdr,fid);
        niftiPadHeader(hdr,fid);
        fclose(fid);
end;

return;


