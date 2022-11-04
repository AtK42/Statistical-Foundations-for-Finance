function ll = tloglik(param, x)
v=param (1); mu=param (2); c=param (3);
if v < 0.01
    v=rand; % An ad hoc way of prevent ing negative values which works , but is NOT recommended! 
end 
if c<0.01
    c=rand;
end
K=beta(v / 2, 0.5) * sqrt(v); z=(x-mu)/c;
ll = -log(c) - log(K) -((v+1)/2) * log(1 + (z.^2) / v);
ll = -sum(ll);