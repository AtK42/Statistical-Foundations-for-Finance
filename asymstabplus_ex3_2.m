function [ f ,F] = asymstabplus_ex3_2 ( xvec , a1 , a2)
% pdf and, optionally, cdf of the asymmetric stable. See also asymstabpdf.m
% Nolan's routine: Call stablepdf(xvec, [a ,b ] ,1)

bordertol=1e-8; lo=bordertol ; hi=1-bordertol ; tol =1e-7;
xl=length(xvec ); F=zeros(xl ,1 ); f=F;
for loop=1:length(xvec )
    %x=(xvec(loop)-d)/c;
    x=xvec(loop);
    f(loop)=(integral(@(uvec) fffplus_ex3_2(uvec,x,a1,a2,1), lo, hi)/pi);
    f(loop)=integral(@(uvec) fffplus_ex3_2(uvec,x,a1,a2,1), lo, hi)/pi;
    if nargout>1
        F(loop)=0.5-(1/pi)*integral(@(uvec) fffplus_ex3_2(uvec,x,a1,a2,0), lo, hi)/pi;
    end
end