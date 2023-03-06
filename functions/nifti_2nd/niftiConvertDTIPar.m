function niftiConvertDTIPar(parName,outName);
%obselete
%niftiConvertDTIPar converts PARREC volume to nifti format
%syntax:
%niftiConvertDTIPar(parName,outName);
%

par = loadPAR_untouch(parName);
%
inPlaneMatrix= par.scn.recon_res;

q = [par.img.special];

bFactor = sort(unique([q.diffusion_b_factor]));
numBFactor = numel(bFactor);
[gradID, gI,gJ] = unique([q.gradient_orientation_number]);
grad = [q(gI).diffusion_ap_fh_lr]';
%sort gradient ID;

numGrad = numel(gradID);
numSlice = par.max.num_slices;

numDirections = numGrad-1;%minus (0,0,0)


meanDWIID = zeros(1,0);
b0ID = zeros(1,0);
otherID = zeros(1,0);
for i = 1:numel(q)
    if(all(q(i).diffusion_ap_fh_lr==0)&&q(i).diffusion_b_factor~=0)
        meanDWIID(end+1) = i;
    elseif(q(i).diffusion_b_factor==0)
        b0ID(end+1) = i;
    else
        otherID(end+1) = i;
    end;
end;

[p,n,e] = fileparts(parName);
recName = fullfile(p,[n,'.rec']);
img = loadREC_ut(recName,par);
%flip y direction
img = flipdim(img,2);


info = [par.img.info];
%
b0Sort = [info(b0ID).slice_num]';
[b0Sorted, b0SortedID ] = sortrows(b0Sort);
b0IDSorted = b0ID(b0SortedID);

b0Image = reshape(img(:,:,1,b0IDSorted),[inPlaneMatrix(1),inPlaneMatrix(2),numSlice]);


if(~isempty(meanDWIID))
    %process Mean Diffusion kurtosis
    DWISort = [toColumnVector([q(meanDWIID).diffusion_b_factor]),toColumnVector([info(meanDWIID).slice_num])];
    [DWISortedd,DWISortedID] = sortrows(DWISort);
    meanDWIIDSorted = meanDWIID(DWISortedID);
    meanDWIImage = reshape(img(:,:,1,meanDWIIDSorted),[inPlaneMatrix(1),inPlaneMatrix(2),numSlice,numBFactor-1]);
end;



if(~isempty(otherID))
    %process true diffusion kurtosis tensor
    %sort images by b values, then gradient directions then slice ID;
    mySort = [[q(otherID).diffusion_b_factor];[q(otherID).gradient_orientation_number];[info(otherID).slice_num]]';
    [mySorted,sortedID] = sortrows(mySort);
    otherIDSorted = otherID(sortedID);
    otherBImage = reshape(img(:,:,1,otherIDSorted),[inPlaneMatrix(1),inPlaneMatrix(2),numSlice,numDirections,numBFactor-1]);

end;

hdr = niftiHeaderFromPar(par);
hdr.dime.dim(1) = 4;
hdr.dime.dim(5) = numDirections*(numBFactor-1)+1;
if(~isempty(meanDWIID))
    hdr.dime.dim(5) = hdr.dime.dim(5)+(numBFactor-1);
end;
hdr = niftiSetImageType(hdr,outName);
niftiSaveHeader(hdr,outName);

fid = niftiOpenImageFileForWrite(outName,hdr);
fid = niftiSaveImageSegment(fid,b0Image,hdr);
if(~isempty(otherID))
    fid = niftiSaveImageSegment(fid,otherBImage,hdr);
end;
if(~isempty(meanDWIID))
    fid = niftiSaveImageSegment(fid,meanDWIImage,hdr);
end;
niftiCloseImageFile(fid);
return;
