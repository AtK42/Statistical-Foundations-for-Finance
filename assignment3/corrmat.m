function B = corrmat(d)
% function to create correlation matrix (positive definite)
% works quite well for low dimensions (until 7 dimensions), otherwise it 
% takes a while
r = -1 + (1+1)*rand(d);
M = tril(r,-1);
B = M + M' + eye(d);

while all(eig(B) > 0) == 0
    r = -1 + (1+1)*rand(d);
    M = tril(r,-1);
    B = M + M' + eye(d);
    disp(1)
end




