function oImg = ifft12c(inImg)
%This function performs centered fft in 1st and 2nd dimension
%function oImg = fft12c(inImg)

oImg = ifftc(ifftc(inImg,[],1),[],2);