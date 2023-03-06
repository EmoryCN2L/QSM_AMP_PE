function X = fftnc(x,siz)

% FFTNC   N-dimensional discrete Fourier Transform with FFT-shift.
%
% X = fftnc(x,siz)
%
% ####### Jian Zhang, jianz@stanford.edu, 08/11/2006 #######

if nargin == 2
    x_s = fftshift(x);
    X_s = fftn(x_s,siz);
    X = fftshift(X_s);
else
    x_s = fftshift(x);
    X_s = fftn(x_s);
    X = fftshift(X_s);
end