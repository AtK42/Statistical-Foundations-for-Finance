function [param, stderr, iters, loglik, Varcov] = GAtestimation(x, vlo, initvec, fixd)
if nargin<2 , vlo = [] ; end, if isempty(vlo), vlo = 0.01; end
if nargin<3 , initvec = []; end, if nargin <4, fixd = []; end
if isempty(initvec)
    versuch=3; vhi=4; loglik = -Inf; vvec=linspace(vlo + 0.02, vhi, versuch);
    for i = 1:versuch
        vv = vvec(i);
        if isempty(fixd), initvec = [2 vv 0.98 0 3];
        else initvec = [vv 0.98 0 3];
        end
        [param0, stderr0, iters0, loglik0, Varcov0] = GAtestimation(x, vlo, initvec, fixd);
        if loglik0 > loglik
            loglik = loglik0; param=param0; stderr=stderr0;
            iters=iters0; Varcov=Varcov0;
        end
    end
    return
end
if isempty(fixd) %  d   v  theta mu  c
    bound.lo =    [0.1 vlo 0.2 -1 1e-4];
    bound.hi =    [ 30 100   3  2 1e+4];
    bound.which = [  1   1   1  0    1];
else %           v theta mu   c
    bound.lo = [vlo 0.2 -1 1e-4];
    bound.hi = [100 3 2 1e+4];
    bound.which = [1 1 0 1];
end

nobs = length(x); maxiter = length(initvec)*100;
tol = 1e-8; MaxFunEvals=length(initvec)*400;
opts=optimset('Display', 'none', 'Maxiter', maxiter, 'TolFun', tol, 'TolX', ...
    tol, 'MaxFunEvals', MaxFunEvals, 'LargeScale', 'Off');
if 1==1
    [pout, fval, exitflag, theoutput, grad, hess] = ...
        fminunc(@(param) GAtloglik(param, x, fixd, bound), ...
            einschrk(initvec, bound), opts);
else
    [pout, fval, exitflag, theoutput] = ...
        fminsearch(@(param) GAtloglik(param, x, fixd, bound), ...
        einschrk(initvec, bound), opts);
    hess=eye(length(pout));
end
V=inv(hess)/nobs; [param, V] = einschrk(pout, bound, V);
param=param'; Varcov=V;
stderr=sqrt(diag(V)); loglik = -fval*nobs;
iters = theoutput.iterations;

function ll = GAtloglik(param, x, fixd, bound)
if nargin<4, bound=0; end
if isstruct(bound), paramvec=einschrk(real(param), bound, 999);
else paramvec=param;
end
if isempty(fixd)
    d=paramvec(1); v=paramvec(2); theta=paramvec(3);
    mu=paramvec(4); c=paramvec(5);
else
    d= fixd; v=paramvec(1); theta=paramvec(2);
    mu=paramvec(3); c=paramvec(4);
end
z=(x-mu)/c; pdf = GAt(z, d, v, theta)/c;
llvec = log(pdf); ll =-mean(llvec);
if isinf(ll), ll =1e5; end