function Dipole_fitting(iname,maskName,oname,opt)
%%This function implements spherical mean filtering for phase imaging
%function SM_filter(iname,maskname, oname,r)
%       iname,      input image name
%       maskname,   mask image name
%       oname,      output image name
%       r,          radius of sphere
%       opt
%           plane, image acquisition plane assume axial
%Deqiang Qiu. qiudeqiang@gmail.com
%Stanford.


%% input
if(~exist('r','var'))
    r = 3;
end
TRUNC = 0.05;
if(~exist('opt','var'))
    opt = [];
end

opt = setoptdefaults(opt,{'writeSF',0,...
                           'B0',3,...
                           'TE',16e-3,...
                           'plane','axial'});
                       
%opt.plane = niftiQFormRot(iname);

%iname = 'P93184_phase_unwrap2.img';
%iname = 'P93184_lapunwrap.img';
%iname = 'test.img';
%maskname = 'P93184_mag_mask_ero5_ero3.img';
%iname = 'phase_unwrap.nii';
%maskname = 'mag_mask_clean5.img';
%oname = 'phase_unwrap_smv.img';

alpha = - 2*pi*42.58 *opt.B0* opt.TE;%this doesn't really matter, but used to control tolerance
%% create convolution kernel
%get image dimension and pixel size
hdr = spm_vol(iname);
orgSz = hdr.dim(1:3);
iRawHdr = niftiLoadHeader(iname);
%pixdim = [hdr.mat(1,1),hdr.mat(2,2),hdr.mat(3,3)];
pixdim = iRawHdr.dime.pixdim(2:4);

%read in images
orginImg = spm_read_vols(hdr)/alpha;
maskImg = spm_read_vols(spm_vol(maskName));
%
paddedSz = ceil(orgSz/2)*2*2;%round image size to even number, then pad the image to double the size

leftpad = floor((paddedSz-orgSz)/2);
rightpad = paddedSz-orgSz-leftpad;

inImg = mypadarray(orginImg,leftpad,rightpad,0);
maskImgPad = mypadarray(maskImg,leftpad,rightpad,0);
%clear('orginImg');


B0_dir = opt.plane(3,:);
if (norm(B0_dir,1)<1.01)
    dipoleFT = fftshift(dipole_kernel( ceil(orgSz/2)*2, pixdim, B0_dir));
else
    dipoleFT = fftshift(dipole_kernel( ceil(orgSz/2)*2, pixdim, B0_dir, 'imagespace'));
end

D2 = fftn(ifftshift(padarray(ifft123c(dipoleFT),(paddedSz-ceil(orgSz/2)*2)/2)));%create non-centric kernel
%D2 = fftn(ifftshift(padarray(ifft123c(GenerateDipoleFT3Drot(ceil(orgSz/2)*2,pixdim,opt.plane)),(paddedSz-ceil(orgSz/2)*2)/2)));%create non-centric kernel


%% conjugate gradient DP fitting
A = FuncMatrixAdaptor(@(x) backwardDF(forwardDF(x,maskImgPad,D2),maskImgPad,D2), @(x) forwardDF(backwardDF(x,maskImgPad,D2),maskImgPad,D2));

b = backwardDF(inImg,maskImgPad,D2);

x = zeros(paddedSz);

num_iter = opt.num_iter;
%figure(100);imagesc(inImg(:,:,round(paddedSz(3)/2)),[-0.2,0.2]);colormap(gray);drawnow
%tic
for i = 1:1
    x = conjugateGradient(A,b,x,1e-10,num_iter,1);

%    dp_fitted = forwardDF(x,maskImgPad,D2);
%    dp_res = inImg-dp_fitted;
%    dp_res = myunpadarray(dp_res,leftpad,rightpad);
%    figure(100);imagesc(dp_res(:,:,round(orgSz(3)/2)),[-0.2,0.2]);colormap(gray);drawnow
end
%toc



    dp_fitted = forwardDF(x,maskImgPad,D2);
    dp_res = inImg-dp_fitted;
    dp_res = myunpadarray(dp_res,leftpad,rightpad);
%    figure(100);imagesc(dp_res(:,:,round(orgSz(3)/2)),[-0.2,0.2]);colormap(gray);drawnow



%% displayimage
%imshowflat(resImg.*maskImg*1000,8,[]);
%imshowflat(dcnvImg*1000,8,[]);

%% output image
[p, n, e] = fileparts(oname);

hdrOut = hdr;
hdrOut.fname = oname;
spm_write_vol(hdrOut,dp_res*alpha);


if(opt.writeSF)
    hdrRes = hdr;
    hdrRes.fname = fullfile(p,[n,'_fld',e]);
    spm_write_vol(hdrRes,myunpadarray(x,leftpad,rightpad));
end

return
