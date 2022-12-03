% define input parameters
% % get random sample of a Student t dist
true_df = 4; dim = 10; n_samp = 20;
rng(42, 'twister');
X = trnd(true_df, dim, n_samp);

% % define the weights of the samples, where each weight must be larger zero and we need sum(w) = 1 (see p. 81)
w = 1/n_samp * ones(n_samp, 1);

% % set the step algorithm, here MMF but you can alternatively define your own function handle
step_algorithm = 'MMF';

% % set maximum number of iterations
anz_steps = 5;

% % set stopping criteria, if 1 then stopping criteria is applied
stop = 1;
abs_criteria = 1;

% % set whether sigma should be regularized to prevent singularity
regularize = 0;

% % set whether negative log-lik should be safed in each step, if yes performence will suffer
save_obj = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mu, nu, sigma, num_steps, time, objective] = iterate_studentT(X, w, step_algorithm, anz_steps, stop, abs_criteria, regularize, save_obj);