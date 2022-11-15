function ES_sim = Simulated_ES_NCT(n_samp, df, mu, alpha, seed)

% function to calculate the ES of a NCT distribution via Simulation. 
% Returns a value for the ES

% params:
% % n_samp: number of samples for the simulation
% % df: degreesof freedom for the NCT
% % mu: Noncentrality parameter for NCT
% % alpha: probability level for ES
% % seed: random seed for reproducability

X = asymtrnd([n_samp 1], mu, df, seed);
VaR = quantile(X, alpha);
Plo = X(X < VaR);
ES_sim = mean(Plo);
end