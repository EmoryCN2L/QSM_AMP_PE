function o = niftiLoadNii(name);
%function o = niftiLoadNii(name);
%output o.hdr,o.img

o.hdr = niftiLoadHeader(name);
o.img = niftiLoadImage(o.hdr,o.hdr.ft,o.hdr.fp,o.hdr.machine);

return;