function img = niftiScaleImage(img,hdr);
%function img = niftiScaleImage(img,hdr);
%scale the image according to slope and intercep specified in the image;

if(hdr.dime.scl_slope==0&&hdr.dime.scl_inter==0)
  img = double(img);
else
  img = double(img)*hdr.dime.scl_slope+hdr.dime.scl_inter;
end

return;
