function oImg = fft123c(inImg)
%This function performs centered fft in 1st and 2nd dimension
%function oImg = fft12c(inImg)

oImg = fftc(fftc(fftc(inImg,[],1),[],2),[],3);