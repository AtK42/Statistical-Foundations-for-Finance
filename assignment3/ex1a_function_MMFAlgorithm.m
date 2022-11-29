function ret = ex1a_function_MMFAlgorithm(true_df, initial_df, reps, dim, n_samp, wgts)
% function to implement the Multivariate Myriad Filter (MMF) as given 
%   on page 91f of "Alternatives to the EM algorithm for ML estimation of
%   location, scatter matrix, and degree of freedom of the Student t
%   distribution (https://doi.org/10.1007/s11075-020-00959-w) by
%   Hasannasab et al.

% input parameters
% % true_df             true degrees of freedom
% % initial_df          starting value for the degrees of freedom
% % reps                number of repetitions
% % dim                 dimension of each random sample
% % n_samp              number of random samples
% % wghts               weights

% output
% % estimate of the degrees of freedom of a Student t distribution

% check input
if initial_df <= 0
    error ("df_initial must be strictly larger 0")
end

if n_samp < dim+1
    error("number of samples must be less than dim + 1 (where dim: sample size)")
end

if sum(wgts == 0) > 0
    error("all weights must be strictly positive")
end

if sum(wgts) ~= 1
    error("the weights must sum to one")
end

% get random sample of a Student t dist
rng(rand(1) * 1000, 'twister');
x_mat = trnd(true_df, dim, n_samp);

% initialization
% % nu (df)
nu_vec = zeros(reps, 1);
nu_vec(1) = initial_df;
% % mu
mu_mat= zeros(dim, reps);
mu_mat(:,1) = sum(x_mat, 2)/n_samp;
% % Sigma
Sigma_mat = zeros(dim, dim, reps);
temp = zeros(dim);
for i = 1:n_samp
    temp = temp + ( x_mat(:,i) - mu_mat(:,1) ) * ( x_mat(:,i) - mu_mat(:,1) )';
end
Sigma_mat(:,:,1) = temp/n_samp;
delta_mat = zeros(n_samp);
gamma_mat = zeros(n_samp);
% looping
for r = 1:reps-1
    % e-step
    for i = 1:n_samp
        delta_mat(i, r) = ( x_mat(:,i) - mu_mat(:,r) )' / Sigma_mat(:,:,r) * ( x_mat(:,i) - mu_mat(:,r) );
        gamma_mat(i, r) = ( nu_vec(r) + dim ) / ( nu_vec(r) + delta_mat(i, r) );
    end
end

ret = nu_vec(end);






