% Function for Question 5, simulated ES for sum of two stable r.v.'s

function ES_sum_sim = Simulated_ES_sum(nobs, a1, a2, b, c, d, xi, seed)
X1 = stabgen(nobs, a1, b, c, d, seed); X2 = stabgen(nobs, a2, b, c, d, seed+1); S = X1 + X2;
q = quantile(S, xi);
ES_sum_sim = [];
for i = 1:3
    Plo = S(S < q(i));
    ES_sum_sim(i) = mean(Plo);
end
end