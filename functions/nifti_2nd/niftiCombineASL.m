function niftiCombineASL(lname,uname,oname,option)
%niftiCombineASL combine ASL stack
%niftiCombineASL(lname,uname,oname,option)
%1name, name of lower stack
%uname, name of upper stack
%oname, output name
%option, specify how to combine the two stack
%        'rl' remove the highest slice from lower stack,
%        'ru',remove the lowest slice from higher stack,
%        'av', average the above two slices

if(exist('option','var')==0)
  option='rl';
end;

if((~strcmpi(option,'av'))&&(~strcmpi(option,'rl'))&&(~strcmpi(option,'ru')))
  error('Invalid option');
end;
    

[lhdr, lft, lfp,lmachine] = load_nii_hdr_ut(lname);
[uhdr, uft, ufp,umachine] = load_nii_hdr_ut(uname);

if (~isequal(lhdr.dime.dim(2:3), uhdr.dime.dim(2:3)))
    
    error('Incompatical inplane matrix: %s dimensions: %s, %s dimension: %s',...
                   lname,num2str(lhdr.dime.dim(2:3)),...
                   uname,num2str(uhdr.dime.dim(2:3)));
    
end;

limg = niftiScaleImage(load_nii_img(lhdr, lft, lfp, lmachine),lhdr);
uimg = niftiScaleImage(load_nii_img(uhdr,uft,ufp,umachine),uhdr);

ohdr = lhdr;
ohdr = niftiSetDataType(ohdr,'float32');
ohdr.dime.scl_inter= 0;
ohdr.dime.scl_slope = 1;
ohdr.dime.dim(4) = lhdr.dime.dim(4)+uhdr.dime.dim(4)-1;

%copy images:

oimg = zeros(ohdr.dime.dim(2),ohdr.dime.dim(3),ohdr.dime.dim(4));
oimg(:,:,1:(lhdr.dime.dim(4)-1)) = limg(:,:,1:(end-1));
oimg(:,:,(lhdr.dime.dim(4)+1):end) = uimg(:,:,2:end);
if(strcmpi(option,'av'))
    oimg(:,:,lhdr.dime.dim(4)) = (limg(:,:,end)+uimg(:,:,1))/2;
elseif(strcmpi(option,'rl'))
    oimg(:,:,lhdr.dime.dim(4)) = uimg(:,:,1);
else %option = 'ru'
    oimg(:,:,lhdr.dime.dim(4)) = limg(:,:,end);
end;

niftiSaveNii(ohdr,oimg,oname,lft);

return;