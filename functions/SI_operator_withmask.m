function A = SI_operator_withmask(kernl,dim,mask)
%SI operator, transform SI image to delta field image, 
%and mask it with mask
%function A = SI_operator_withmask(kernl,dim,mask)
%kernl, dipole kernel, e.g. generated using GenerateDipoleFT3D
%dim, image dimension
%mask, mask image
%
%Deqiang Qiu. qiudeqiang@gmail.com
%c)Stanford. 2010

    A.times = @(x) SI_op_times(x,kernl,dim,mask);
    A.trans = @(x) SI_op_trans(x,kernl,dim,mask);
return


function y = SI_op_times(x,kernl,dim,mask)
    num_vox = prod(dim);
    y = ifftnc(fftnc(reshape(x,dim)).*kernl);
    if(isempty(mask))
        y = y(:);
    else
        y = y(mask>0);
    end
    
return

function y = SI_op_trans(x,kernl,dim,mask)
    num_vox = prod(dim);
    tmpX = zeros(dim);
    if(isempty(mask))
        tmpX = x;
    else
        tmpX(mask>0) = x;
    end
    y = reshape(ifftnc(fftnc(tmpX).*kernl),[num_vox,1]);
    return