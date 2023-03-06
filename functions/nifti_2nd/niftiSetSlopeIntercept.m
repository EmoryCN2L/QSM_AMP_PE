function hdr = niftiSetSlopeIntercept(hdr,scl)
%Set slope and intercept of a nifti Header

%
hdr.dime.scl_slope = scl(1);
hdr.dime.scl_inter = scl(2);

return;