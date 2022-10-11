% Function for Question 5, simulated ES for sum of two stable r.v.'s

function ES_sum_sim = Simulated_ES_sum(nobs, a1, a2, b, c, d, xi, seed)
nobs = 1e6;
X1 = stabgen(nobs, a1, b, d, c, seed); X2 = stabgen(nobs, a2, b, d, c, seed); S = X1 + X2;
q = quantile(S, xi);
ES_sum_sim = [];
for i = 1:3
    Plo = S(S < q(i));
    ES_sum_sim(i) = mean(Plo);
end
end