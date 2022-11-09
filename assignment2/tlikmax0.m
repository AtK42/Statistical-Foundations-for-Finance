function MLE = tlikmax0(x, initvec)
% Program Listing 4.5: Attempts to maximize the log-likelihood of the i.i.d. Student’s t
% model usingdata vector x, starting values initvec=[df location scale], aconvergence
% tolerance of 0.00001 and allowing at most 200 function evaluations.
% x: data
% initvec = [df location scale]
tol =1e-5;
opts=optimset('Disp', 'none', 'LargeScale', 'Off', ...
    'TolFun', tol, 'TolX', tol, 'Maxiter',200);
MLE = fminunc(@(param) tloglik(param, x), initvec, opts);
    
function ll = tloglik(param, x)
v=param (1); mu=param (2); c=param (3);
if v < 0.01
    v = rand; % An ad hoc way of preventing negative values which works , but is NOT recommended! 
end 
if c < 0.01
    c = rand;
end
K = beta(v / 2, 0.5) * sqrt(v); z=(x-mu)/c;
ll = -log(c) - log(K) -((v+1)/2) * log(1 + (z.^2) / v);
ll = -sum(ll);