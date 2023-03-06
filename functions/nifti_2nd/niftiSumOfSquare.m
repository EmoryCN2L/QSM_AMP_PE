function niftiSumOfSquare(inames,oname)
%

hdr = spm_vol(inames{1});
ohdr = hdr;
oimg = zeros([hdr.dim]);

for i = 1:numel(inames)
    img = spm_read_vols(spm_vol(inames{i}));
    oimg = oimg + img.^2;
end
oimg = sqrt(oimg);
ohdr.fname = oname;

spm_write_vol(ohdr,oimg);
