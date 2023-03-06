function mat = niftiQFormRot(hdr)
%This function calculates the sform transformation matrix 

if(ischar(hdr))
    hdr = niftiLoadHeader(hdr);
end

if(hdr.dime.pixdim(1)>=0)
    qfac = 1;
else
    qfac = -1;
end

ss = sqrt(sum([hdr.hist.quatern_b^2,hdr.hist.quatern_c^2,hdr.hist.quatern_d^2]));
if(ss>1 && ss<1+1e-4);%forgive rounding error
    hdr.hist.quatern_b = hdr.hist.quatern_b/ss;
    hdr.hist.quatern_c = hdr.hist.quatern_c/ss;
    hdr.hist.quatern_d = hdr.hist.quatern_b/ss;
end
qt = eye([4,4]);
qt(1:3,1:3) = quaterion2mat(hdr.hist.quatern_b,hdr.hist.quatern_c,hdr.hist.quatern_d);

qt(1:3,4) = [hdr.hist.qoffset_x;hdr.hist.qoffset_y;hdr.hist.qoffset_z];


mat = qt*diag([1,1,qfac,1]);
mat = mat(1:3,1:3);

return;