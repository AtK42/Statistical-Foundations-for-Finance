function ES_sim = Simulated_ES_symStable(n_samp, tail_index, scale, loc, alpha, seed)

% function to calculate the ES of a (symmetric) stable distribution via Simulation. 
% Returns a value for the ES

% params:
% % n_samp: number of samples for the simulation
% % tail_index: tail index, conventionally denoted by alpha
% % alpha: probability level for ES
% % seed: random seed for reproducability

X = stabgen(n_samp, tail_index, 0, scale, loc, seed);
VaR = quantile(X, alpha);
Plo = X(X < VaR);
ES_sim = mean(Plo);
end