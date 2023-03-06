function outData = diffusionModelKurtosis(par, bvals)

lnS = par(1);
D = par(2);
K = par(3);
outData = lnS-bvals*D+1/6*bvals.^2*D.^2*K;

return;