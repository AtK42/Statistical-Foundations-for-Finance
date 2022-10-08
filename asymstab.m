function [ f ,F] = asymstab ( xvec , a , b )
% pdf and, optionally, cdf of the asymmetric stable. See also asymstabpdf.m
% Nolan's routine: Call stablepdf(xvec, [a ,b ] ,1)

if nargin<3, b=0; end
bordertol=1e-8; lo=bordertol ; hi=1-bordertol ; tol =1e-7;
xl=length(xvec ); F=zeros(xl ,1 ); f=F;
for loop=1:length(xvec )
    x=xvec(loop);
    f(loop)=integral(@(uvec) fff(uvec,x,a,b,1), lo, hi)/pi;
    if nargout>1
        F(loop)=0.5-(1/pi)*integral(@(uvec) fff(uvec,x,a,b,0), lo, hi)/pi;
    end
end