function v=vech(M,invert)
if nargin>1 % M is a vector, and we return a matrix
    vec=M; n=(sqrt(8*length(vec)+1)-1)/2;
    V=zeros(n,n);
    for i=1:n
        take=(n-i+1); comp=vec(1:take); vec=vec(take+1:end);
        V(i:n,i)=comp; V(i,i:n)=comp';
    end, v=V;
else % M is a matrix, and we return vech(M)
    tt=tril(M); v=tt(tt~=0); %tril: Lower triangular part of matrix
end