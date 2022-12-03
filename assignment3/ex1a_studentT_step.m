% This code belongs to the paper
%
% M. Hasannasab, J. Hertrich, F. Laus, and G. Steidl. 
% Alternatives to the EM algorithm for ML-estimation of location, scatter
% matrix and degree of freedom of the student-t distribution.
% Numerical Algorithms, 2020.
% DOI: https://doi.org/10.1007/s11075-020-00959-w
%
% If you use this code, please cite the paper.
%
% This function performs one step of the MMF for
% estimating the parameters of a Student-t distribution.
%
% INPUTS:
%     X                 - d x n array containing the samples
%     w                 - 1 x n array containing the weights of the samples
%     nu_r,mu_r,sigma_r - current estimates of nu, mu and sigma
%     regularize        - we add in each step regularize*eye(d) to sigma to
%                         prevent that sigma becomes singular.
%                         If X contains at least n affine independent samples,
%                         this parameter should be set to 0.
%     nu_step           - declares the nu-step of the algorithm.
%                         Choose nu_step=@nu_step_mmf for the MMF.
%     delta_r           - 1 x n array containing the values of delta_r. Set
%                         delta_r=nan to compute the values within this
%                         function.
%
% OUTPUTS:
%     nu_r_plus_one     - updated estimate of nu.
%     mu_r_plus_one     - updated estimate of mu.
%     sigma_r_plus_one  - updated estimate of sigma.
%     delta_r_plus_one  - updated values of delta.
function [nu_r_plus_one, mu_r_plus_one, sigma_r_plus_one, delta_r_plus_one]=ex1a_studentT_step(X, w, nu_r, mu_r, sigma_r, regularize, nu_step, delta_r)

[d,n]=size(X);

% set default value
if nargin<7 || any(isnan(delta_r(:)))
    delta_r=sum(((sigma_r^(-1))*(X-repmat(mu_r,1,n))).*(X-repmat(mu_r,1,n)),1);
end

% M-step: update the parameters
% % for mu
mu_r_plus_one=sum(repmat(w.*(1./(nu_r+delta_r)),d,1).*X, 2)/sum(w.*(1./(nu_r+delta_r)), 2);
%mu_r_plus_one=sum(repmat(w'.*(1./(nu_r+delta_r)),d,1).*X, 2)/sum(w'.*(1./(nu_r+delta_r)), 2);
% % for sigma
sigma_r_plus_one=((X-repmat(mu_r_plus_one,1,n)).*repmat(w./(nu_r+delta_r),d,1))*(X-repmat(mu_r_plus_one,1,n))'/sum(w .*(1./(nu_r+delta_r)),2)+regularize*eye(d);
%sigma_r_plus_one=((X-repmat(mu_r_plus_one,1,n)).*repmat(w'./(nu_r+delta_r),d,1))*(X-repmat(mu_r_plus_one,1,n))'/sum(w'.*(1./(nu_r+delta_r)),2)+regularize*eye(d);
% % for delta
delta_r_plus_one=sum(((sigma_r_plus_one^(-1))*(X-repmat(mu_r_plus_one,1,n))).*(X-repmat(mu_r_plus_one,1,n)),1);
[nu_r_plus_one,~]=nu_step(nu_r,delta_r_plus_one,d,n);

% fail save constraint
if nu_r_plus_one>1e+5
   nu_r_plus_one=1e+5; 
end

end