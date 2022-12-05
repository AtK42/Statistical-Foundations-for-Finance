function [final_nu, nu_vec, mu, sigma, stop] = ex1a_function_MMFAlgorithm_ver3(x_mat, initial_df, wgts, reps, progress_update)
% function to implement the Multivariate Myriad Filter (MMF) as given 
%   on page 91f of "Alternatives to the EM algorithm for ML estimation of
%   location, scatter matrix, and degree of freedom of the Student t
%   distribution (https://doi.org/10.1007/s11075-020-00959-w) by
%   Hasannasab et al.
%
% input parameters
% % x_mat               matrix of random samples of a multivarite t dist
% % initial_df          starting value for the degrees of freedom
% % wgts                weights
% % reps                number of repetitions
%
% output
% % nu_vec              estimate of the degrees of freedom of a Student t distribution
% % mu                  mean vector of the latest iteration
% % sigma               variance-covariance matrix of the latest iteration


% check input
if initial_df <= 0
    error ("df_initial must be strictly larger 0")
end

if size(x_mat, 1) > size(x_mat, 2)
    error("number of samples must be less than dim + 1 (where dim: sample size)")
end

if sum(wgts == 0) > 0
    error("all weights must be strictly positive")
end

if sum(wgts) - 1 > 1e-10
    error("the weights must sum to one")
end

% start timing
tic

% initialize variables
% % get dimension of the data
d = size(x_mat , 1);
% % nu
nu = initial_df;
nu_vec = zeros(reps, 1);
% % mu
%mu = 1/length(x_mat)*sum(x_mat, 2);
mu = sum(x_mat, 2)/size(x_mat, 2);
%initialize sigma matrix
sigma0 = 0;
for i = 1:length(x_mat)
    sigma0 = sigma0 + (x_mat(:, i) - mu)*(x_mat(:, i) - mu)';
end
sigma = sigma0/length(x_mat);

if progress_update == 1
    disp('********************************');
end
for r = 1:reps
    % E-step: compute weights
    % % initialize vectorss for delta and gamma
    delta = zeros(size(x_mat, 2), 1);
    gamma = zeros(size(x_mat, 2), 1);
    % % fill vectors with values for the current rep loop
    for i = 1:length(x_mat)
       delta(i) = (x_mat(:, i) - mu)' / sigma * (x_mat(:, i) - mu);
       gamma(i) = (nu + d)/(nu + delta(i,1));
    end

    % M-step: update the parameters
    % % initialize (set to zero) variable to save denominator for updating mu and sigma
    denom = 0;
    for i = 1:size(x_mat, 2)
        denom = denom + wgts(i)*gamma(i);
    end

    % % mu
    % % % initialize (set to zero) variable to save the numerator
    mu_nom = 0;
    % % % calculate the nominator for updating mu
    for i = 1:size(x_mat, 2)
        mu_nom = mu_nom + wgts(i) * gamma(i) * x_mat(:, i);
    end
    % % % update mu value
    mu_lagged = mu;
    mu = mu_nom / denom;

    % % sigma
    % % % initialize (set to zero) variable to save the numerator
    sigma_num = 0;
    % % % calculate the nominator for updating sigma
    for i = 1:length(x_mat)
        sigma_num = sigma_num + ( wgts(i) * gamma(i) * (x_mat(:,i) - mu) * (x_mat(:,i) - mu)' );
    end
    % % % update sigma value
    sigma_lagged = sigma;
    sigma = sigma_num / denom;

    % % nu
    % % % initialize (set to zero) variable to save the sum part of the updating step for nu
    nu_sum = 0;
    for i = 1:length(x_mat)
        nu_sum = nu_sum + wgts(i) * ( (nu + d) / (nu + delta(i)) - log( (nu + d)/(nu + delta(i)) ) - 1 );
    end
    nu = fzero(@(x)  phi_func(x/2) - phi_func((x+d) / 2) + nu_sum , [1e-100, 1e100]);
    nu_vec(r) = nu;
    
    % convergence criteria
    nu_stop = false; mu_stop = false(size(x_mat, 1), 1); sigma_stop = false(size(x_mat, 1));
    if r > 2 && nu_vec(r) - nu_vec(r-1) < 1e-10 && nu_vec(r-1) - nu_vec(r-2) < 1e-10
        nu_stop = true;
        for i = 1:size(x_mat, 1)
            if mu(i) - mu_lagged(i) < 1e-10
                mu_stop(i) = true;
            end
            for j = 1:size(x_mat, 1)
                if sigma(i,j) - sigma_lagged(i,j) < 1e-10
                    sigma_stop(i,j) = true;
                end
            end
        end
    end
    if nu_stop == true && sum(mu_stop) == size(x_mat, 1) && sum(sum(sigma_stop)) == (size(x_mat, 1))^2
        stop = r;
        break
    end

    % user friendly feature for progress updates
    if progress_update == 1
        if mod(r, 50) == 0
            disp(['finished ', num2str(r), ' reps']);
        end
    end
end % r-loop (reps)
final_nu = nu;

% end timing
time = toc;

if progress_update == 1
    disp('*********** RESULTS  ***********');
    disp(['__for sample size = ', num2str(size(x_mat, 2)), '__'])
    disp('estimated df: '); disp(final_nu);
    disp('estimated mu vector: '); disp(mu);
    disp('estimated sigma matrix: '); disp(sigma);
    disp(['time to run: ', num2str(time), 's']);
    disp(['iterations completed: ', num2str(stop)]);
    disp('********************************')
end
end % function

function [phi] = phi_func(x)
    % see p. 81
    phi = psi(x) - log(x);
end