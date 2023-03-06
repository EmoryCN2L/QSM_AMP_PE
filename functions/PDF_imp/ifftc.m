function o = ifftc(varargin)
%This function performs centered fft
%function o = ifftc(img,N,dim)
%
%Deqiang Qiu, Stanford University

img = varargin{1};
if(nargin>=2)
    N = varargin{2};
else
    N = [];
end
if(nargin>=3)
    dim = varargin{3};
else
    dim = [];
end

if(isempty(dim))
    o = ifftshift(ifftn(fftshift(img)));
else
    o = ifftshift(ifft(fftshift(img,dim),N,dim),dim);
end
