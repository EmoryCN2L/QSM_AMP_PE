function o = GenerateDipoleFT3Drot(dim,pixsize,plane)
%This function generates 3D dipole kernal based on the dimension and pixel
%size of the 3 dimensions as well as matrix plane
%
%function o = generateDipoleFT3D(dim,pixsize,plane)
%dim, 3-element vector of the dimension
%pixsize, 3-element pixel size
%plane, image plane, 'axial'(default),'cor','sag',or 3x3 rotation matrix
%
%Deqiang Qiu. qiudeqiang@gmail.com
%c)Stanford. 2011

o = zeros(dim);

dk = 1./(dim.*pixsize);%k space element size

org = (dim)/2+1;%considering fftshift, the origin should be at (dim)/2+1;

if(~exist('plane','var'))
    plane = 'axial';
end

rot = zeros(3,3);
if(ischar(plane)) %only the last row matters
    if(strcmpi(plane,'axial'))
        %zidx = 3;
        rot(3, 3) =1;
    elseif(strcmpi(plane,'cor'))
        %zidx = 2;
        rot(3, 2) = 1;
    elseif(strcmpi(plane,'sag'))
        rot(3, 2) = 1; 
        %zidx = 1;
    else
        error('unrecognized plane');
    end
elseif(ismatrix(plane))
    if(~isequal(size(plane),[3,3]))
        error('rotation matrix has to be 3 by 3');
    end
    
    %if(abs(abs(det(plane))-1)>1e-3)
    %    error('rotation matrix should det=1 or -1');
    %end
    rot = plane;
else
    error('unrecognized plane/rotation matrix');
end

zproj = rot(3,:);


for idx1 = 1:dim(1);
    for idx2 = 1:dim(2);
        for idx3 = 1:dim(3);
            idx = [idx1,idx2,idx3];
            if(all(idx==org))
                o(idx1,idx2,idx3)  =0;
            else
                kz = sum(zproj.*((idx-org).*dk));
                o(idx1,idx2,idx3) = 1/3-kz^2/sum(((idx-org).*dk).^2);
            end
        end
    end;
end


return



