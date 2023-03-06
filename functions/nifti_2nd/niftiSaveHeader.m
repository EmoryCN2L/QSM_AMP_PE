function hdr = niftiSaveHeader(hdr,fileName,filetype)
%function niftiSaveHeader(hdr,fileName,filetype)
%hdr,
%fileName,
%filetype, optional
%Deqiang Qiu.
%
[prefix,type] = niftiImagePrefix(fileName);

if(~isfield(hdr,'machine'))
    hdr.machine = 'ieee-le';
end

if(exist('filetype','var')>0)
    type=filetype;
end;

switch type;
    case 0,
        hdr.hist.magic = '';
        hdr.dime.vox_offset = 0;
        hdrName = [prefix,'.hdr'];
        fid = fopen(hdrName,'wb',hdr.machine);
        assert(fid~=-1,sprintf('Fail to open file %s for writing',hdrName));
        niftiSaveHeaderFid(hdr,fid);
        fclose(fid);
    case 1,
        hdr.hist.magic = 'ni1';
        hdr.dime.vox_offset = 0;
        hdrName = [prefix,'.hdr'];
        fid = fopen(hdrName,'wb',hdr.machine);
        assert(fid~=-1,sprintf('Fail to open file %s for writing',hdrName));
        niftiSaveHeaderFid(hdr,fid);
        fclose(fid);
    case 2,
        hdr.hist.magic = 'n+1';
        if(hdr.dime.vox_offset<352)
            hdr.dime.vox_offset = 352;
        end
        hdrName = [prefix,'.nii'];
        fid = fopen(hdrName,'wb',hdr.machine);
        assert(fid~=-1,sprintf('Fail to open file %s for writing',hdrName));
        fseek(fid,0,'bof');
        niftiSaveHeaderFid(hdr,fid);
        niftiPadHeader(hdr,fid);
        fclose(fid);
end;

return;


