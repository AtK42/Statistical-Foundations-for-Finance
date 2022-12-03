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
% This function implements the nu-step of the MMF.
% INPUTS:
%       nu_r            - current nu
%       delta_r         - 1 x n array, where delta_r(i)=(x_i-mu_r)'*sigma_r^(-1)*(x_i-mu_r)
%       dim_dat         - dimension of the data
%       n_samp          - number of samples
% OUTPUTS:
%       nu_r_plus_one   - next nu
%       evals           - number of evaluations of the subfunctions A and B
function [nu_r_plus_one,evals]=ex1a_nu_step_mmf(nu_r, delta_r, dim_dat, n_samp)

% define function (first part, \phi(...) on p. 92, Algo 3) for calculation of \nu_{r+1} using digamma funciton and log (see p. 81)
A = @(x)psi(x/2)-log(x/2);

% define derivative of function A (see line above, resp. p. 92 for A)
derivative_A = @(x).5*psi(1,x/2)-1./x;

% define the sum (second part) for calculation of \nu_{r+1} on p. 92
B = @(x)sum((x+dim_dat)./(x+delta_r) - log((x+dim_dat)./(x+delta_r)) -1)/n_samp;

% apply above defined function on the current instance of nu_r
b_nu_r = B(nu_r);

% define entire function, i.e., put the two pieces (A and B) from above
% together and additionally its derivative
% % note that we need the derivative for Newton Raphson below (which uses
% % x_{n+1} = x_n - f(x_n)/derivative_f(x_n), n \in \mathbb{N})
f = @(nu)A(nu)-A(nu+dim_dat)+b_nu_r;
derivative_f = @(nu)derivative_A(nu)-derivative_A(nu+dim_dat);

% use Newton-Raphson for determining the zeroes in the nu-steps
[nu_r_plus_one, evals] = ex1a_newton(nu_r, f, derivative_f);

end