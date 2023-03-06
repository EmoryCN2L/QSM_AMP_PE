function hdr = niftiHeaderFromPar(par)
%create nifti header from par file
%set fields: image matrix and slice resolution
%assume 1 for the 4 the dimension
%
 %%%TODO calculate real rotation information from angulation

q = [par.img.orient];
orient = unique([q.slice_orientation]);
if(numel(orient)~=1)
    error('Multiple slice orietation in the image');
end;
if(orient~=1)
    warning('Non transverse orientation not well tested');
end;
    numSlice = par.max.num_slices;
    hdr = niftiHdr();
    [sMat,realVoxelSize,quaternXYZ,quaternBCD,qFac] = parToNiftiMat(par.orient.ang_midslice,par.orient.off_ctr_midslice,[par.scn.recon_res,par.max.num_slices],par.scn.fov,orient);
    hdr.dime.dim(1) = 4;
    hdr.dime.dim(2:5) = [par.scn.recon_res(1),par.scn.recon_res(1),numSlice,1];
    hdr.dime.pixdim(1) = qFac;
    %to do check acqusition plane, assume axial
    hdr.dime.pixdim(2:4) = realVoxelSize;
    %print a warning for non axial acquistion
    hdr.hist.qform_code = 1;
    hdr.hist.sform_code = 1;
    hdr.hist.srow_x = sMat(1,:);
    hdr.hist.srow_y = sMat(2,:);
    hdr.hist.srow_z = sMat(3,:);
    hdr.hist.quatern_b = quaternBCD(1);
    hdr.hist.quatern_c = quaternBCD(2);
    hdr.hist.quatern_d = quaternBCD(3);
    hdr.hist.qoffset_x = quaternXYZ(1);
    hdr.hist.qoffset_y = quaternXYZ(2);
    hdr.hist.qoffset_z = quaternXYZ(3);
    
return;