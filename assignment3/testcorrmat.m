counting = zeros(10,1);
for i = 1:10
    tic
    [a, b] = corrmat(8);
    toc
    counting(i) = b;
end
mean(counting)
