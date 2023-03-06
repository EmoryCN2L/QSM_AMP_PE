function fid = niftiOpenImageFileForWrite(fileName,hdr,filetype)
%function fid = niftiOpenImageFileForWrite(fileName,hdr,filetype)
%
%if 
[prefix, type] = niftiImagePrefix(fileName);
if(exist('filetype','var')>0&&~isempty(filetype))
    type=filetype;
end;

switch type;
    case {0,1},
        imgName= [prefix,'.img'];
        fid = fopen(imgName,'wb',hdr.machine);
        assert(fid~=-1,sprintf('Fail to open file %s for writing',imgName));
    case 2
        imgName= [prefix,'.nii'];
        fid = fopen(imgName,'r+b',hdr.machine);
        assert(fid~=-1,sprintf('Fail to open file %s for writing',imgName));
        
        %do neccessary padding
        %move to the end of the file
        if(hdr.dime.vox_offset<352)
            hdr.dime.vox_offset =352;
        end
        
        fseek(fid,0,'eof');
        pos = ftell(fid);
        if(pos<hdr.dime.vox_offset)
            niftiPadHeader(hdr,fid);
        else
            fseek(fid,hdr.dime.vox_offset,'bof');
            %sanity check
            pos = ftell(fid);
            if(pos~=hdr.dime.vox_offset)
                error('Fail to pad the header');
            end;     
        end
            
        
       
end;

