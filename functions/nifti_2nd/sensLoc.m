function locs = sensLoc(loc,mtx, sensFac)
%function locs = sensLoc(mtx, sensFac)
%loc, starts from 1 in the accelerated matrix
%locs, starts from 1 in the extendedMatrix

%this is the reconstructed matrix size
extendSz = ceil(mtx*sensFac);

addedSz = extendSz-mtx;

leftHalf = floor(addedSz/2);
rightHalf = addedSz-leftHalf;

%Center of the matrix
ct = ceil(extendSz/2);

%initialize the array to the non-aliased one
locs = loc+leftHalf;%add the current one

%search to right
idx = loc + leftHalf + mtx;
while(idx<=extendSz)
    locs = [locs, idx];
    idx = idx + mtx;
end

%search to left
idx = loc + leftHalf - mtx;
while(idx>=1)
    locs = [locs, idx];
    idx = idx - mtx;
end

