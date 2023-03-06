function niftiMean(fileNames,oname)
%Calculate the mean value over a series of images. 
%function niftiMean(fileNames,oname)
%The filename can be cell array, or a struct array with "name" field
%


if(isstruct(fileNames)&&isfield(fileNames,'name'));
    inames = {fileNames.name};
else
    inames = fileNames;
end;
consistentImage = 1;
numImages = numel(inames);
isSameSlopeIntercept = 1;


for i = 1:numImages;
    header(i).hdr = niftiLoadHeader(inames{i});
    %check consistency
    if(i ==1)
        numDim = header(1).hdr.dime.dim(1);
        errorMsg = sprintf('Inconsistent Image Dimension:\n Image 1, Dimensions: %s,Pixel spacing:%s \n',...
                            num2str(header(1).hdr.dime.dim(1:(numDim+1))),num2str(header(1).hdr.dime.dim(1:(numDim+1)))); 
    else
        if(~isequal(header(i).hdr.dime.dim(1:(numDim+1)),header(1).hdr.dime.dim(1:(numDim+1)))||...
           ~isequal(header(i).hdr.dime.pixdim(1:(numDim+1)),header(1).hdr.dime.pixdim(1:(numDim+1))))
            consistentImage = 0;
            errorMsg = [errorMsg,sprintf('Image %d: Dimensions: %s , pixel spacing: %s\n',i,...
                                        num2str(header(i).hdr.dime.dim(1:(numDim+1))),num2str(header(i).hdr.dime.dim(1:(numDim+1))))]; 
        end
    end;
end;

if(consistentImage ==0)
    error(errorMsg);
end;

%go on

oHeader = header(1);
oHeader.img = zeros(oHeader.hdr.dime.dim(2:(numDim+1)));

oHeader.hdr.dime.scl_inter = 0;
oHeader.hdr.dime.scl_slope = 1;

for i = 1:numImages;
    oHeader.img = oHeader.img+niftiScaleImage(load_nii_img(header(i).hdr,header(i).hdr.ft,header(i).hdr.fp,header(i).hdr.machine),header(i).hdr)/numImages;
    disp(['working on ',int2str(i)])
end;
minVal = min(oHeader.img(:));
maxVal = max(oHeader.img(:));
minValDataType = niftiType(oHeader.hdr.dime.datatype,'minval');
maxValDataType = niftiType(oHeader.hdr.dime.datatype,'maxval');


if(minVal<minValDataType||maxVal>maxValDataType)
    warning('Untested feature, scaling');
    if((maxVal-minVal)<(maxValDataType-minValDataType))
        oHeader.hdr.dime.scl_slope = 1;
    else
        oHeader.hdr.dime.scl_slope = (maxVal-minVal)/(maxValDataType-minValDataType);
    end;
    oHeader.hdr.dime.scl_inter = minVal-oHeader.hdr.dime.scl_slope*minValDataType;
    oHeader.img = (oHeader.img-minVal)/oHeader.hdr.dime.scl_slope+minValDataType;
end;



oHeader.hdr = niftiSaveHeader(oHeader.hdr,oname);
fid = niftiOpenImageFileForWrite(oname,oHeader.hdr);
niftiSaveImageSegment(fid,oHeader.img,oHeader.hdr);
fclose(fid);

