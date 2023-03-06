function niftiEigenToTensorTrackvis(root,oname,varargin);
%generate tensor image from eigen vectors and eigen values images
%function niftiEigenToTensor(root,oname);
%
%Option: Value pairs
%'postfix', followed by 'fsl', or a cell array of postfix;
%
%
%This output order of tensor element is compatible with trackvis; which is
%strange: D11,D12,D22,D13,D23,D33
%

postfix={'_v1.nii','_v2.nii','_v3.nii','_e1.nii','_e2.nii','_e3.nii'};

i = 1;
while i<=numel(varargin)
  if(strcmpi(varargin{i},'postfix'))
    i = i+1;
    if(isstr(varargin{i})&&strcmpi(varargin{i},'fsl'))
      postfix={'_V1.nii','_V2.nii','_V3.nii','_L1.nii','_L2.nii','_L3.nii'};
    elseif(iscell(varargin{i}))
      postfix=varargin{i};
    else
      error('Unrecognized postfix format');
    end;
  end;
  
  i = i+1;
  
end;



%check all file exist
for i= 1:numel(postfix);
  assert(exist([root,postfix{i}],'file')~=0,'%s not existi',[root,postfix{i}]);
end;

img = cell(numel(postfix),1);

%convert eigen vectors/values to tensor
for i = 1:numel(postfix);
  [hdr(i),ft{i},fp{i},machine{i}] = load_nii_hdr_ut([root,postfix{i}]);
  img{i} = niftiScaleImage(load_nii_img(hdr(i),ft{i},fp{i},machine{i}),hdr(i));
end;


tensor = eigenToTensor(img{1},img{2},img{3},img{4},img{5},img{6},'trackvis');
%reorder tensor

hdrOut = hdr(1);

hdrOut.dime.scl_slope = 1;
hdrOut.dime.scl_inter = 0;
hdrOut.dime.dim(1) = 5;
hdrOut.dime.dim(5) = 1;
hdrOut.dime.dim(6) = 6;

%change storage type to single float
hdrOut.dime.datatype = 16;
hdrOut.dime.bitpix=32;
%chage to .nii filetype
hdrOut.hist.magic= 'n+1';
hdrOut.dime.vox_offset = 352;


fidOut = fopen(oname,'wb',machine{1});
assert(fidOut~=-1,'Fail to open %s for writing',oname);

save_nii_hdr(hdrOut,fidOut);

%padHeader
curPos =ftell(fidOut);
if(curPos<hdrOut.dime.vox_offset)
  bytewrote = fwrite(fidOut,zeros(hdrOut.dime.vox_offset-curPos,1),'uchar');%pad header
end;


write_nii_segment(fidOut,tensor,hdrOut);


return;