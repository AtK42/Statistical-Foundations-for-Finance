function[param,stderr,iters,loglik,Varcov]=tlikmax(x,initvec)

%%%%%%%%      df    mu    c
bound.lo=   [  .1     0     0.01  ];
bound.hi=   [100   100   100     ];
bound.which=[  1     1     1     ];
%In this case , as bound . which for mu is zero , mu will not be
%restricted. As such, the values for .lo and .hi are irrelevant
%listing 4.6

maxiter = 100; tol = 1e-3;% change these as you see fit
opts = optimset('Display', 'notify-detailed' , 'Maxiter', maxiter, 'TolFun', tol, 'TolX', tol, 'LargeScale', 'Off');
[pout, fval, exitflag, theoutput, grad, hess] = fminunc(@(param) tloglik(param, x, bound), einschrk(initvec, bound), opts);

V = inv(hess); % Don't negate: we work with the neg of the loglik
[param, V] = einschrk(pout, bound, V); % Transform back, apply delta method
param = param' ; Varcov = V;
stderr = sqrt(diag(V)); % Approx std err of the params
loglik = -fval; % The value of the loglik at its maximum.
iters = theoutput.iterations; % Number of loglik function evals

function ll=tloglik(param, x ,bound)
if nargin < 3
    bound=0;
end
if isstruct(bound)
    paramvec=einschrk(real(param), bound, 999);
else
    paramvec = param;
end

v = paramvec(1); mu = paramvec(2);c = paramvec(3);
K = beta(v/2, 0.5) * sqrt(v); z = (x - mu)/c;
ll = -sum(-log(c) - log(K) - ((v+1)/2) * log(1+ (z.^2)/v));