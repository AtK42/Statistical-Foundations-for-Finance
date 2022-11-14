function[df, ncp, loc, scale] = nctlikmax(x, initvec, pdfopt)
% MLE of the singly NCT distribution for all four parameters, being
% the degrees of freedom (df), the noncentrality parameter (ncp), the 
% location (loc) and the scale (scale)

% as input, the data (x) is needed and an initial vector as the starting
% point for the estimation (initvec = [df ncp loc scale]))

% additional parameter to determine what method to use to calculate the pdf
% according to the exercise: 1 is the one from matlab and everything else
% is his approximation method

tol =1e-5; maxiter = 1000;
opts=optimset('Disp', 'none', 'LargeScale', 'Off', ...
    'TolFun', tol, 'TolX', tol, 'Maxiter',200);
%opts = optimset('Display', 'notify-detailed' , 'Maxiter', maxiter, 'TolFun', tol, 'TolX', tol, 'LargeScale', 'Off');
%opts = optimoptions('fminunc', 'Display', 'notify-detailed' , 'MaxIterations', maxiter, 'OptimalityTolerance', tol, 'StepTolerance', tol);
MLE = fminunc(@(param) nctloglik(param, x, pdfopt), initvec, opts)
    
function ll = nctloglik(param, x, pdfopt)
df=param(1); ncp=param(2); loc=param(3); scale=param(4);
if df < 0.01 df = rand; end 
if scale < 0.01 scale = rand; end
z=(x-loc)/scale;
if pdfopt == 1
    pdf = nctpdf(z, df, ncp)/scale; %from matlab
else
    pdf = exp(stdnctpdfln_j(z, df,ncp))/scale; %his approximation
end
llvec = log(pdf); ll=-mean(llvec);
if isinf(ll), ll =1e5; end
