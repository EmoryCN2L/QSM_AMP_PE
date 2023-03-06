function nifti2DICOMByTemplate(niftiFile,oDir,hdr)
%% convert nifti file to DICOM images and use a dicom template file for other 
%function nifti2DICOMByTemplate(niftiFile,oDir,dcmTemplate)
%
%Deqiang Qiu, qiudeqiang@gmail.com, 2015

%hdr = dicominfo(dcmTemplate);
im = niftiLoadNii(niftiFile);


%% get origin and x, y, z vector from nift field
RAS2LPS = diag([-1,-1,1]);
rotMat = niftiQFormRot(im.hdr)*diag(im.hdr.dime.pixdim(2:4));
org = RAS2LPS*[im.hdr.hist.qoffset_x; im.hdr.hist.qoffset_y; im.hdr.hist.qoffset_z];
xvec = RAS2LPS*rotMat(1:3,1);
yvec = RAS2LPS*rotMat(1:3,2);
zvec = RAS2LPS*rotMat(1:3,3);

rotMat4x4 = [xvec,yvec,zvec,org];
rotMat4x4(4,1:4) = [0,0,0,1];

%normalize them
norm_xvec = xvec./sqrt(xvec'*xvec);
norm_yvec = yvec./sqrt(yvec'*yvec);
norm_zvec = zvec./sqrt(zvec'*zvec);
%% retrieve the scale factor
intercept = im.hdr.dime.scl_inter;
slope = im.hdr.dime.scl_slope;
if(slope==0)
    slope=1;
    intercept = 0;
end

img = im.img*slope + intercept;

%%
MINVAL = 0;
MAXVAL = 4095;
imgMin = min(img(:));
imgMax = max(img(:));
% 
% inter = imgMin;
% slp = (imgMax-imgMin)/MAXVAL;
% 
% img_rescale = (img-inter)/slp;
%%

seriesUid = dicomuid;
%%loop and output dicom images

for idx = 1:im.hdr.dime.dim(4)
    oName = fullfile(oDir,sprintf('IM-%d.dcm',idx));
    headerDCM = hdr;
    headerDCM.Rows = im.hdr.dime.dim(3);
    headerDCM.Columns = im.hdr.dime.dim(2);
    position = rotMat4x4*[0;0;idx-1;1];
    sliceLocation = norm_zvec'*position(1:3);
    
    headerDCM.ImagePositionPatient = position(1:3);
    headerDCM.ImageOrientationPatient = [norm_xvec;norm_yvec];
    headerDCM.SliceThickness = im.hdr.dime.pixdim(4);
    headerDCM.PixelSpacing = im.hdr.dime.pixdim(2:3)';
    headerDCM.SliceLocation = sliceLocation;
    
    %instance number
    headerDCM.InstanceNumber = idx;
    headerDCM.SeriesInstanceUID = seriesUid;
    headerDCM.SOPInstanceUID = dicomuid;
    headerDCM.MediaStorageSOPInstanceUID = headerDCM.SOPInstanceUID;
    headerDCM.WindowCenter = round(imgMin/2+imgMax/2);
    headerDCM.WindowWidth = round(imgMax-imgMin);
%      headerDCM.RescaleSlope = slp;
%      headerDCM.RescaleIntercept = inter;
%      headerDCM.RescaleType = 'US';
    dicomwrite(int16(img(:,:,idx)'),oName,headerDCM);
end

%%
return;




