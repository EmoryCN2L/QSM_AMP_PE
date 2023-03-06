function niftiSincInterpOneDim(iname,oname,dim,interpLength);
%
%

[hdr, ft, fp, machine ] = load_nii_hdr_ut(iname);
if(interpLength == hdr.dime.dim(1+dim))
    %TODO copy files
    return;
end;

img = load_nii_img(hdr,ft,fp,machine);
orgSize = size(img);
numDim = numel(orgSize);

%put the dimension of interest to the first dimension
img = shiftdim(img,dim-1);
s = size(img);
shiftSize = zeros(numel(s),1);
shiftSize(1) = round(s(1)/2);

fftImg = circshift(fft(img,[],1),shiftSize);

offtImg = zeros([interpLength,s(2:end)]);
leftPadSize = floor((interpLength-s(1))/2);
offtImg((leftPadSize+1):(leftPadSize+s(1)),:,:,:) = fftImg;

backShift = zeros(numel(s),1);
backShift(1) = -1*(shiftSize(1)+leftPadSize);

oImg = shiftdim(ifft(circshift(offtImg,backShift),[],1),numDim-dim+1);

ohdr = hdr;
ohdr.dime.dim(1+dim) = interpLength;
ohdr.dime.pixdim(1+dim) = hdr.dime.pixdim(1+dim)*hdr.dime.dim(1+dim)/ohdr.dime.dim(1+dim);

fid = fopen(oname,'wb');
save_nii_hdr(ohdr,fid);
niftiPadHeader(ohdr,fid);
write_nii_segment(fid,oImg,ohdr);
fclose(fid);

