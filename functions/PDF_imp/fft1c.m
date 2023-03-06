function oImg = fft1c(inImg)
%This function performs centered fft in 1st and 2nd dimension
%function oImg = fft12c(inImg)

oImg = fftc(inImg,[],1);
