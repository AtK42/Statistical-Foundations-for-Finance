% program listing 12.2
function y = mvtpdfmine(x,df,mu,Sigma)
% x is a d X 1 vector. Unlike Matlab's version, cannot pass a matrix.
% Matlab's routine accepts a correlation (not dispersion) matrix.
% So, just need to do the usual scale transform. For example:
% x=[0.2 0.3]'; C = [1.4; .4 1]; df = 2;
% scalevec=[1 2]'; xx=x./scalevec; mvtpdf(xx,C,df)/prod(scalevec)
% Same as:
% Sigma = diag(scalevec) * C * diag(scalevec); mvtpdfmine(x,df,[],Sigma)
d=length(x);
if nargin<3, mu = []; end, if isempty(mu), mu = zeros(d,1); end
if nargin<4, Sigma = eye(d); end
x = reshape(x,d,1); mu = reshape(mu,d,1); term = (x-mu)' * inv(Sigma) * (x-mu);
logN=-((df+d)/2)*log(1+term/df); logD=0.5*log(det(Sigma))+(d/2)*log(df*pi);
y = exp(gammaln((df+d)/2) - gammaln(df/2) + logN - logD);