function [x,r] = conjugateGradient(A,b,x0,tol,maxit,printFreq)
%This function finds the solution to Ax = b through conjugate gradient method
%function [x,r] = conjugateGradient(A,b,x0,tol,maxit,printFreq)
%   x0,     initial value
%   tol,    iteration stops when relative residue is below this tolerance level
%   maxit,  maximum number of iterations
%   printFreq, controls how often we print out info
%
%By defaults this function draw an example of iteration paths when the
%problem is two dimensional for education purposes, since probably you
%don't need to use CG for this case. You can turn this off by setting
%variable "illustration" to 0.
%
%Deqiang Qiu, qiudeqiang@gmail.com
%Stanford University, 2011

if(~exist('tol','var'))
    tol = 1e-6;
end
if(~exist('maxit','var'))
    maxit = 1000;
end

if(~exist('x0','var'))
    x0 = b*0;
end

if(~exist('printFreq','var'))
    printFreq = 50;
end


%%
dim = size(b,1);
maxit = min(dim,maxit);

illustration = 1;

if(dim>2)
    illustration = 0;
end;


%% 2x2 illustration

if(illustration)
    [X, Y] = meshgrid([-6:0.1:6],[-4:0.1:4]);
    Z = X*0;
    for i = 1:numel(X)
        Z(i) = 1/2*[X(i),Y(i)]*A*[X(i);Y(i)]-b'*[X(i);Y(i)];
    end
    figure;
    title('Conjugate Gradient');
    contour(X,Y,Z,100);
    hold on;
end

%%
r = b-A*x0;
x = x0;
d = r;

res = tol+1;
it = 1;

while(res>tol && it <=maxit)
    x_ = x;
    sigma_old = r(:)'*r(:);
    
    Ad = A*d;
    
    a = sigma_old/(d(:)'*Ad(:));
    
    x = x+a*d;
    if(mod(it,50)==0)%recalcualte residual to correct for rounding errors
        r = b-A*x;
    else
        r = r-a*Ad;
    end
    sigma_new = r(:)'*r(:);
    beta = sigma_new/sigma_old;
    d = r+beta*d;%form conjugate 
    
    
    res = norm(x(:)-x_(:),2)/norm(x(:),2);
    %do 2D plog
    if(illustration)
        plot([x_(1),x(1)],[x_(2),x(2)],'o-');
    end
    
    if(mod(it,printFreq)==0)
        fprintf('Interation %d, residual%f \n',it,sigma_new);
    end
    it = it+1;
end

return;