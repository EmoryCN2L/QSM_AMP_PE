

%%
hdr = niftiLoadHeader('template.img');
ohdr = hdr;


img = zeros([hdr.dime.dim(2:4), 3,62]);

for iec = 1:5
    for v = 1:62
        for is = 1:3
            for sl = 1:17
                im = dicomread(sprintf('e08403s012-v%03d-shot%02d-echo%02d-%02d.mag.dcm',v,is,iec,sl));
                img(:,:,sl,is,v) = flipdim(permute(im,[2,1]),2);
            end
        end
    end
    ohdr.fname = sprintf('magecho%d.img',iec);
    ohdr.dime.dim(1) = 4;
    ohdr.dime.dim(5) = 3*62;
    oImg = reshape(img,ohdr.dime.dim(2:5));
    niftiSaveNii(ohdr,oImg,ohdr.fname);
    
end