function f = convopdf(svec, a1, a2)
%function that integrates the convointegrand function at different values
%of s (given as input: svec) with respect to x
xl=length(svec); f=zeros(xl, 1);
for loop=1:length(svec)
    s=svec(loop);
    f(loop)=integral(@(xvec) convointegrand(xvec,s,a1,a2), -Inf, Inf);
end
