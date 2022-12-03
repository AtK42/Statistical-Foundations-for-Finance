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
% This function implements the algorithm MMF for Student-t distributions 
% using Newton-Raphson to find the roots.
%
% INPUTS:
%     X             - d x n array containing the samples
%     w             - 1 x n array containing the weights of the samples
%     step_algorithm- Choice of the algorithm, which determines the nu, mu
%                     and sigma step. step_algorithm has to be of type
%                     string.
%                     For step_algorithm='MMF' the function performs the
%                     MMF.
%                     For step_algorithm='GMMF' the function performs the
%                     GMMF.
%                     For step_algorithm='aEM' the function performs the
%                     accelerated EM algorith,
%                     For step_algorithm='ECME' the function perfors the
%                     ECME algorithm.
%                     Alternatively the step_algorithm can be declared by a
%                     function handle.
%     anz_steps     - Maximum number of iterations
%     stop          - stop==1 means we apply the stopping criteria, 
%                     stop==0 means we perform exactly anz_steps spteps.
%     abs_criteria  - Replaces the relative stopping criteria by an
%                     absolute stopping criteria
%     regularize    - we add in each step regularize*eye(d) to sigma to
%                     prevent that sigma becomes singular.
%                     If X contains at least n affine independent samples,
%                     this parameter should be set to 0.
%     save_obj      - save_obj==1 means we save in each step the negative
%                     log-likelihood value. Otherwise save_obj should be
%                     set to 0. Note that save_obj==1 slows down the
%                     performance of the algorithm.
%
% OUTPUTS:
%     mu,nu,sigma   - estimated parameters
%     num_steps     - if stop==1 then num_steps contains the number of
%                     steps, until the stopping criteria is reached.
%                     Otherwise NaN is returned.
%     time          - execution time of the algorithm
%     objective     - if save_obj==1, objective(r) contains the negative
%                     log-likelihood value of nu_r,mu_r,Sigma_r.
%                     Otherwise objective(r)=0.
function [mu, nu, sigma, num_steps, time, objective]=ex1a_iterate_studentT(X, w, step_algorithm, anz_steps, stop, abs_criteria, regularize, save_obj)

% set some default values
if nargin<8
    save_obj=0;
end
if nargin<7
    regularize=0;
end
if nargin<6
    abs_criteria=0;
end
if nargin<5
    stop=0;
end

% check and define function handle
% % note that on the line where step_algorithm is being defined, the 
% % function ex1a_nu_step_mmf is needed
if ~isa(step_algorithm,'function_handle')
    if isa(step_algorithm,'char')
        step_algorithm=@(X,w,nu_r,mu_r,sigma_r,regularize,delta_r)ex1a_studentT_step(X, w, nu_r, mu_r, sigma_r, regularize, @ex1a_nu_step_mmf, delta_r);
    else
        error('Unknown algorithm. step_algorithm has to be a function handle or a string.')
    end
end
tol=1e-5;

% set the tolerance if absolute stopping criteria has been chosen
if abs_criteria==1
    tol=1e-8;
end

% initialization
[dim,n_samp]=size(X);
L=@(nu,mu,sigma,w,delta)-2*gammaln((dim+nu)/2)+2*gammaln(nu/2)-nu.*log(nu)+(dim+nu).*sum(w.*log(nu+delta))+log(det(sigma));
objective=zeros(1,anz_steps+1);
mu_r=sum(X,2)/n_samp;
sigma_r=(X-repmat(mu_r,1,n_samp))*(X-repmat(mu_r,1,n_samp))'/n_samp+regularize*eye(dim);
nu_r=3; % set to arbitrary value

% start timing
tic;

% E-step
arg_cont=1;
r=0;
delta_r=sum(((sigma_r^(-1))*(X-repmat(mu_r,1,n_samp))).*(X-repmat(mu_r,1,n_samp)),1);
while arg_cont==1
    % counter
    r=r+1;
    
    % check for numerical errors, i.e. negative delta_r values and NAs
    if sum(delta_r<0)>0 || isnan(nu_r)
        error('Numerical errors')
    end
    
    % save log-lik if so chosen
    if save_obj==1
        if any(isnan(delta_r))
            delta_r=sum(((sigma_r^(-1))*(X-repmat(mu_r,1,n_samp))).*(X-repmat(mu_r,1,n_samp)),1);
        end
        objective(r)=L(nu_r, mu_r, sigma_r, w, delta_r);
    end
    
    % use ex1a_studenT_step to get the next iterations of nu, mu, sigma and delta
    [nu_r_plus_one, mu_r_plus_one, sigma_r_plus_one, delta_r_plus_one]=step_algorithm(X, w, nu_r, mu_r, sigma_r, regularize, delta_r);
    
    % check whether algo should be stopped or continue
    if stop==1 && abs_criteria==0
        eps=((sum((mu_r-mu_r_plus_one).^2)+sum(sum((sigma_r-sigma_r_plus_one).^2)))/(sum(mu_r.^2)+sum(sum(sigma_r.^2))))^.5+((log(nu_r_plus_one)-log(nu_r))^2/log(nu_r)^2)^.5;
        if eps<tol
            arg_cont=0;
        end
        if r==10000
            arg_cont=0;
        end
    end
    
    if stop==1 && abs_criteria==1
        eps=(sum((mu_r-mu_r_plus_one).^2)+sum(sum((sigma_r-sigma_r_plus_one).^2))+(log(nu_r_plus_one)-log(nu_r))^2)^.5;
        if r==10000
            arg_cont=0;
        end
        if eps<tol
            arg_cont=0;
        end
    end
    
    % reassign values for mu, sigma, nu and delta to be the current values
    mu_r=mu_r_plus_one;
    sigma_r=sigma_r_plus_one;
    nu_r=nu_r_plus_one;
    delta_r=delta_r_plus_one;
    if stop==0
        if r==anz_steps
            arg_cont=0;
        end
    end
end

% save log-lik if desired
%if save_obj==1
%    delta_here=sum(((sigma_r^(-1))*(X-repmat(mu_r,1,n_samp))).*(X-repmat(mu_r,1,n_samp)),1);
%    objective(r+1)=L(nu_r,mu_r,sigma_r,w,delta_here);
%end

% stop the clock
time=toc;
if stop==1
    num_steps=r;
else
    num_steps=nan;
end

% return values for mu, nu and sigma
mu=mu_r; nu=nu_r;sigma=sigma_r;
end