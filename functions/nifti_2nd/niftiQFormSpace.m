function mat = niftiQFormSpace(hdr)
%This function calculates the sform transformation matrix 

if(ischar(hdr))
    hdr = niftiLoadHeader(hdr);
end

if(hdr.dime.pixdim(1)>=0)
    qfac = 1;
else
    qfac = -1;
end

qt = eye([4,4]);
qt(1:3,1:3) = quaterion2mat(hdr.hist.quatern_b,hdr.hist.quatern_c,hdr.hist.quatern_d);

qt(1:3,4) = [hdr.hist.qoffset_x;hdr.hist.qoffset_y;hdr.hist.qoffset_z];


mat = qt*diag([hdr.dime.pixdim(2),hdr.dime.pixdim(3),qfac*hdr.dime.pixdim(4),1]);

return;