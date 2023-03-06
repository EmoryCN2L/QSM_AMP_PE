function mat = niftiSFormRot(hdr)
%calculate SForm from a nifti header
%

%offset = [1, 0, 0, hdr.hist.qoffset_x;
%          0, 1, 0, hdr.hist.qoffset_y;
%          0, 0, 1, hdr.hist.qoffset_z;
%          0, 0, 0, 1];

qRot = niftiQuaternion2Mat([hdr.hist.quatern_b, hdr.hist.quatern_c, hdr.hist.quatern_d]);
qFac = (hdr.dime.pixdim(1)>=0)*1+(hdr.dime.pixdim(1)<0)*-1;
qFacMat = diag([1,1,qFac,1]);
%voxMat = diag([hdr.dime.pixdim(2:4),1]);

mat = qRot*qFacMat;
mat = mat(1:3,1:3);
return;
