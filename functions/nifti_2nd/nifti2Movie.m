function nifti2Movie(iname, oname,fps,lim)
%This function converts nifti file or 3D matrix to movie
%function nifi2Movie(iname, oname)
%Syntax:
%
%Input:
%   iname, input file name or matrix
%   oname, output file name
%   fps, frame per second
%   lim, image window level
%Output
%
%Deqiang Qiu. qiudeqiang@gmail.com
%Stanford University, 2011
if(ischar(iname))
    img = spm_read_vols(spm_vol(iname));  
else
    img = iname;
end

hf = figure;
ha = gca;
if(~exist('lim','var'))
    lim = [];
end

if(isempty(lim))
    lim = [min(img(:)),max(img(:))];
end
    

himg = image(squeeze(img(:,:,1)'),'parent',ha,'cdatamapping','scaled');
set(gca,'clim',lim);
set(gca,'ydir','normal');
colormap gray;
set(gca,'xtick',[],'ytick',[]);

aviobj = avifile(oname);
for i = 1:size(img,3)
    set(himg,'CData',squeeze(img(:,:,i)'));
    frm(i) = getframe(hf);
    aviobj=  addframe(aviobj,frm(i));
end

close(hf);
aviobj = close(aviobj);

return


