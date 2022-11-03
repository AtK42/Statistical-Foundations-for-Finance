function [f,F] = asymstab_book(xvec,a,b)
    bordertol=1e-8; lo=bordertol; hi=1-bordertol; tol=1e-6;
    xl=length(xvec); F=zeros(xl,1); f=F;
    for loop=1:length(xvec)
        x=xvec(loop); dopdf=1;
        f(loop)= quadl(@fff,lo,hi,tol,[],x,a,b,1) / pi;
        if nargout>1
            F(loop)=0.5-(1/pi)* quadl(@fff,lo,hi,tol,[],x,a,b,0);
        end
    end
end