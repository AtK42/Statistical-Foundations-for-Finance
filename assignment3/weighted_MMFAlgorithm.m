function [final_nu, nu_vec, mu, sigma] = weighted_MMFAlgorithm(rho, x_mat, initial_df, reps)
% function to implement the Multivariate Myriad Filter (MMF) combined with
% hyperbolic weight decay as of page 590 in the book Linear Models and time
% series analysis by Marc S. Paolella.
%
% input parameters
% % rho                 parameter for weight decay
% % x_mat               matrix of random samples of a multivarite t dist
% % initial_df          starting value for the degrees of freedom
% % reps                number of repetitions
%
% output
% % nu_vec              estimate of the degrees of freedom of a Student t distribution
% % mu                  mean vector of the latest iteration
% % sigma               variance-covariance matrix of the latest iteration

% calculate weight vector
T = length(x_mat); tvec = (1:T); omega=(T-tvec + 1).^(rho-1); w = omega'/sum(omega);
disp(sum(w));

% call MMF algorithm
[final_nu, nu_vec, mu, sigma] = ex1a_function_MMFAlgorithm_ver3(x_mat, initial_df, w, reps)
end