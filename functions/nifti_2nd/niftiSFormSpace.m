function mat = niftiSFormSpace(hdr);
%This function calculates the sform transformation matrix 

if(ischar(hdr))
    hdr = niftiLoadHeader(hdr);
end

mat = eye([4,4]);
mat(1,:) = hdr.hist.srow_x;
mat(2,:) = hdr.hist.srow_y;
mat(3,:) = hdr.hist.srow_z;

return