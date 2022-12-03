%% 1a
% Implement the "algorithm 3 MMF" in the attached paper. 
% Part of this assignment is to train you to look at research articles and 
%   implement the methods.
% In this case, you do NOT have to read the theory, but rather just code 
%   that algorithm, for which they give very nice pseudo-code. 

% input parameters
%true_df = 4; % freely assumed
%initial_df = 1; % freely assumed
%reps = 10; % freely assumed
%dim = 10; % sample size, freely assumed
%n_samp  = 20; % number of samples, n > d, freely assumed
%wgts = 1/n_samp * ones(n_samp, 1); % each weight must be larger zero and we need sum(wgts) = 1 (see p. 81)

% call function
%[nu_final, delta_mat, gamma_mat, mu_mat, Sigma_mat, nu_vec, x_mat] = ex1a_function_MMFAlgorithm(true_df, initial_df, reps, dim, n_samp, wgts);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% call function
[mu, nu, sigma, num_steps, time, objective] = ex1a_iterate_studentT(X, w, step_algorithm, anz_steps, stop, abs_criteria, regularize, save_obj);
%% 1b
% Simulate a 3-variate IID multivariate Student t (with, say, 4 df), zero 
%   mean vector but please a NON-DIAGONAL Sigma matrix that you invent ---
%   obviously, it has to be positive definite.
% So, another part of this assignment involves you being able to generate 
%   a covariance matrix.
 
% NOTE: I choose a zero location vector because the quality of the 
%   estimation results are probably mostly "invariant" to this choice. 
% As an example, with the linear regression model, the choice of the true 
%   beta regression coefficients you would use in simulation exercises is 
%   irrelevant, in the sense that certain properties of the OLS estimator 
%   and the residuals are invariant to its actual value. (This is because 
%   of the nature of the projection matrix used, but that is a topic for 
%   next semester in regression...)

% Task 1b comment (i).
% You take the main diagonal to be ones, so that you are in fact making a 
%   correlation matrix. But in the below simulations, you do NOT assume it 
%   is known that these values are unity. We just do this so we have nice 
%   numbers, and it is very possible that the performance of the MLE is 
%   (at least closely) invariant to scaling, similar to my assumption on 
%   the location terms.
 
% Task 1b comment (ii).
% For this low dimension of 3, you can trivially generate a valid 
%   correlation matrix with some trial and error. Or, google it, there are 
%   surely routines out there for general construction. Or, for fun 
%   (and maybe extra points), invent your own algorithm.

% for further comments see mail


df = 4;
mean_vec = zeros(3, 1);
%% 1c
% Simulate your 3-d MVT with say T=200 and T=2000 observations, and 
%   estimate it using the "MMF" algorithm. Repeat this 500 times, and 
%   report wonderful colored boxplots, such as the 2nd graphic seen here,
%   https://ch.mathworks.com/help/stats/boxplot.html
%   for your parameter estimates. Thus, we assess the quality of the MLE 
%   (via their algorithm).
 
% Crucially, you SUBTRACT THE TRUE VALUE of the parameter from your 
%   estimates before you boxplot them, so the boxplot shows the DEVIATION 
%   FROM THE TRUTH. Got it?
%% 1d
% Using the SAME 500 replications as above (so results are even more 
%   comparable), do the same but using your own custom made MLE routine 
%   using brute force log likelihood maximization. Notice this is super 
%   easy, since I did it all for you, with 2 dimensions. You just need to 
%   make a simple extension to 3 dimensions (3-d). See Program Listings 
%   12.1 and 12.2 in my time series book. 
 
% More ambitiously, and for some bonus points, and very instructive, make 
%   the 3-d MVT MLE program for dimension general d. This is actually also 
%   easy, because I did it already: See section 12.5.3 in my time series 
%   book for how to do this with matlab for a related model.
 
% How many parameters do you have? 
%   3 for the location vector, 
%   one for df, and 
%   6 for the correlation matrix. 
%   (We assume we do NOT know that the diagonals are unity, so we estimate 
%   them.) So, this is easy for fminunc or fmincon in Matlab for this small 
%   number of parameters.
 
% Be sure to constrain the parameters with at least obvious box constraints 
%   as I do with "einschraenk" or whatever I called it. More advanced is to 
%   use fmincon and have the nonlinear constraint that the Sigma matrix has 
%   to be positive definite.
%% 1e
% Using the SAME 500 replications as above (so results are even more 
%   comparable), do the same but using the ECME algorithm from
%   C Liu and D B Rubin, (1995) "ML estimation of the t distribution using 
%   EM and its extensions, ECM and ECME", Statistica Sinica, [5, pp19-39]
%   http://www3.stat.sinica.edu.tw/statistica/oldpdf/A5n12.pdf
%   and also the fast *approximate* method from:
%   C Aeschlimna, J Park and KA Cak, "A Novel Parameter Estimation Algorithm 
%   for the Multivariate t-Distribution and its Application to Computer 
%   Vision" [ECCV 2010]
%   http://link.springer.com/chapter/10.1007%2F978-3-642-15552-9_43

% YOU DON'T have to program those methods yourself. They are freely 
%   available on the web, for Matlab:
%   https://github.com/robince/tdistfit

% report the usual boxplots also for these methods, so we can see the 
%   accuracy of these methods. We expect that brute force MLE and the EM 
%   algorithm (ECME in this case) and the MMF algorithm from the paper I 
%   sent you, all result in the same parameter estimation results.
% The approximation method from Aeschlimna, Park and Cak presumably will 
%   not be as accurate, but it should be super fast.

% also set up the tic and toc (or the more advanced ways Matlab allows for 
%   timings) and compare and report estimation times. These are only fair 
%   if the requested accuracy is the same across the algorithms, so see if 
%   this is easy to do. 
% So, bottom line: Also report the computation times for the various 
%   algorithms. We obviously expect the brute force MLE to be the slowest.