function [ret, delta_mat, gamma_mat, mu_mat, Sigma_mat, nu_vec, x_mat] = ex1a_function_MMFAlgorithm(true_df, initial_df, reps, dim, n_samp, wgts)
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

if sum(wgts) - 1 > 1e-10
    error("the weights must sum to one")
end

% get random sample of a Student t dist
rng(42, 'twister');
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
delta_mat = zeros(n_samp, reps);
gamma_mat = zeros(n_samp, reps);

% looping
% % e-step
for r = 1:reps-1
    for i = 1:n_samp
        delta_mat(i, r) = ( x_mat(:,i) - mu_mat(:,r) )' / Sigma_mat(:,:,r) * ( x_mat(:,i) - mu_mat(:,r) );
        gamma_mat(i, r) = ( nu_vec(r) + dim ) / ( nu_vec(r) + delta_mat(i, r) );
    end
%end

% % m-step
%for r = 1:reps-1
    temp_mu_num = zeros(dim, 1);
    temp_denom = 0;
    temp_Sigma = zeros(dim, dim);
    temp_sum = 0;
    % % % mu
    for i = 1:n_samp
        temp_mu_num = temp_mu_num + wgts(i) * gamma_mat(i, r) * x_mat(:,i);
        temp_denom = temp_denom + wgts(i) * gamma_mat(i, r);
    end
    mu_mat(:, r+1) = temp_mu_num / temp_denom;
    % % % Sigma
    for i = 1:n_samp
        temp_Sigma = temp_Sigma + wgts(i) * gamma_mat(i, r) * ( x_mat(:, i) - mu_mat(:, r+1) ) * ( x_mat(:, i) - mu_mat(:, r+1) )';
    end
    Sigma_mat(:,:,r+1) = temp_Sigma / temp_denom;
    % % % nu
    for i = 1:n_samp
        temp_sum = temp_sum + wgts(i) * ( (nu_vec(r) + dim) / (nu_vec(r) + delta_mat(i, r+1)) ) - log( (nu_vec(r) + dim) / (nu_vec(r) + delta_mat(i, r+1)) ) -1;
    end

    syms nu_solve
    res = solve( [(psi(nu_solve / 2)-log(nu_solve/2)) - (psi((nu_solve + dim) / 2) - log((nu_solve + dim) / 2)) + temp_sum == 0], nu_solve );

    % only to get some results, delete everything except for line
    % nu_vec(r+1) = ... when things run correctly
    if isempty(res)
        nu_vec(r+1) = 1;
    else
        disp(res);
        nu_vec(r+1) = res; % since solution is unique (o.w. try subs() )
    end
end

ret = nu_vec(end);
%delta_mat = delta_mat;
%gamma_mat = gamma_mat;
%mu_mat = mu_mat;
%Sigma_mat = Sigma_mat;
%nu_vec = nu_vec;
