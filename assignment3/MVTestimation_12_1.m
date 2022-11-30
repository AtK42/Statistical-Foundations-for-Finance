% program listing 12.1
function [param,stderr,iters,loglik,Varcov] = MVTestimation(x)
% param: (k, mu1, mu2, Sigma_11, Sigma_12, Sigma_22)
[nobs d] = size(x); if d~=2, error('not done yet, use EM'), end
if d == 2
    %%%%%%%% k mu1 mu2 s11 s12 s22
    bound.lo = [ 0.2 -1 -1 0.01 -90 0.01];
    bound.hi = [ 20 1 1 90 90 90];
    bound.which = [ 1 0 0 1 1 1];
    initvec = [2 -0.8 -0.2 20 2 10];
end
maxiter=300; tol=1e-7; MaxFunEvals=length(initvec)*maxiter;
opts=optimset('Display','iter','Maxiter',maxiter,'TolFun',tol,'TolX',tol,...
'MaxFunEvals',MaxFunEvals,'LargeScale','Off');
[pout,fval,~,theoutput,~,hess]= fminunc(@(param) MVTloglik(param,x,bound),einschrk(initvec,bound),opts);
V=inv(hess)/nobs; % Don't negate because we work with the negative of the loglik
[param, V] = einschrk(pout, bound, V); % transform and apply delta method to get V
param=param'; Varcov=V; stderr=sqrt(diag(V)); % Approximate standard errors
loglik = -fval*nobs; iters=theoutput.iterations;

function ll = MVTloglik(param,x,bound)
if nargin<3, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
[nobs d]=size(x);
%%% two dimensional case (original see book)
%Sig=zeros(d,d); k=param(1); mu=param(2:3); % Assume d=2
%Sig(1,1)=param(4); Sig(2,2)=param(6); Sig(1,2)=param(5); Sig(2,1)=Sig(1,2);
k=param(1); mu=param(d+2:2*d+1); Rterms=param(3*d+2:end);
Rt=Rterms; RR=[1 , Rt(1:d-1)]; Rt=Rt(d:end);
for i = 2:d-1, RR = [RR, 1, Rt(1:d-i)]; Rt = Rt(d-i+1:end); end
RR=[RR, 1]; R=vech(RR,1);
if min(eig(R)) < 1e-10, ll=1e5;
else pdf=zeros(nobs,1);
for i=1:nobs
    pdf(i) = mvtpdfmine(x(i,:),k,mu,R);
end
llvec=log(pdf); ll=-mean(llvec); if isinf(ll), ll=1e5; end
end