function [B, counter] = corrmat(d)
% function to create correlation matrix (positive definite)
% works quite well for low dimensions (until 7 dimensions), otherwise it 
% takes a while
r = -1 + 2*rand(d);
M = tril(r,-1);
B = M + M' + eye(d);
counter = 1;

while all(eig(B) > 0) == 0
    r = -1 + (1+1)*rand(d);
    M = tril(r,-1);
    B = M + M' + eye(d);
    counter = counter + 1;
end





