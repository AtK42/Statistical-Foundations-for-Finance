function[df, ncp, loc, scale] = nctlikmax(x, initvec)
% MLE of the singly NCT distribution for all four parameters, being
% the degrees of freedom (df), the noncentrality parameter (ncp), the 
% location (loc) and the scale (scale)

% as input, the data (x) is needed and an initial vector as the starting
% point for the estimation (initvec = [df ncp loc scale]))

tol =1e-5; maxiter = 1000;
opts=optimset('Disp', 'none', 'LargeScale', 'Off', ...
    'TolFun', tol, 'TolX', tol, 'Maxiter',200);
%opts = optimset('Display', 'notify-detailed' , 'Maxiter', maxiter, 'TolFun', tol, 'TolX', tol, 'LargeScale', 'Off');
%opts = optimoptions('fminunc', 'Display', 'notify-detailed' , 'MaxIterations', maxiter, 'OptimalityTolerance', tol, 'StepTolerance', tol);
MLE = fminunc(@(param) nctloglik(param, x), initvec, opts)
    
function ll = nctloglik(param, x)
df=param(1); ncp=param(2); loc=param(3); scale=param(4);
if df < 0.01 df = rand; end 
if scale < 0.01 scale = rand; end
z=(x-loc)/scale;
pdf = nctpdf(z, df, ncp)/scale; %from matlab
%pdf = exp(stdnctpdfln_j(z, df,ncp)); %his approximation
llvec = log(pdf); ll=-mean(llvec);
