%% 1a
% Implement the "algorithm 3 MMF" in the attached paper. 
% Part of this assignment is to train you to look at research articles and 
%   implement the methods.
% In this case, you do NOT have to read the theory, but rather just code 
%   that algorithm, for which they give very nice pseudo-code. 

% input parameters
true_df = 4; % freely assumed
initial_df = 5; % freely assumed
iter = 500; % freely assumed
dim = 3; % sample size, freely assumed
%n_samp_low = 2e3; % number of samples, n > d, freely assumed
n_samp_high = 2e4;

%wgts_low = 1/n_samp_low * ones(n_samp_low, 1); % each weight must be larger zero and we need sum(wgts) = 1 (see p. 81)
wgts_high = 1/n_samp_high * ones(n_samp_high, 1); % each weight must be larger zero and we need sum(wgts) = 1 (see p. 81)

% get random sample of a Student t dist
%rng(4, 'twister');
vcov_mat = corrmat(dim);
%x_mat_low = mvtrnd(vcov_mat, true_df, n_samp_low)';
x_mat_high = mvtrnd(vcov_mat, true_df, n_samp_high)';

% call function
%[final_nu_low, nu_vec_low, mu_low, sigma_low, stop_low] = ex1a_function_MMFAlgorithm_ver3(x_mat_low, initial_df, wgts_low, iter, 1);
[final_nu_high, nu_vec_high, mu_high, sigma_high, stop_high] = ex1a_function_MMFAlgorithm_ver3(x_mat_high, initial_df, wgts_high, iter, 1);

% plot
figure
%plot(1:stop_low, nu_vec_low(1:stop_low), '--', 'Color', 'blue', 'DisplayName', ['samp size = ', num2str(n_samp_low)]);
%hold on
plot(1:stop_high, nu_vec_high(1:stop_high), 'Color', 'green', 'DisplayName', ['samp size = ', num2str(n_samp_high)]);
yline(true_df, 'Color', 'red', 'DisplayName', ['true df = ', num2str(true_df)]);
title('convergence of the estimated degrees of freedom using MMF algorithm');
xlabel('iterations');
ylabel('df');
ylim([0, 7]);
%legend(n_samp_low, n_samp_high, 'true df');
legend('Location', 'northeast');
hold off
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
reps = 500; true_df = 4; n_samp_low = 2e2; n_samp_high = 2e3; dim = 3;

% storage
% % variance-covariance matrix
vcov_mat = zeros(dim, dim, reps);
% % random sample
x_mat_low_store = zeros(dim, n_samp_low, reps); x_mat_high_store = zeros(dim, n_samp_high, reps);

% see function corrmat(d)
for r = 1:reps
    % get random sample of a Student t dist
    %rng(4, 'twister');
    vcov_mat(:,:,r) = corrmat(dim);
    x_mat_low_store(:,:,r) = mvtrnd(vcov_mat(:,:,r), true_df, n_samp_low)';
    x_mat_high_store(:,:,r) = mvtrnd(vcov_mat(:,:,r), true_df, n_samp_high)';
end


% 1c
% Simulate your 3-d MVT with say T=200 and T=2000 observations, and 
%   estimate it using the "MMF" algorithm. Repeat this 500 times, and 
%   report wonderful colored boxplots, such as the 2nd graphic seen here,
%   https://ch.mathworks.com/help/stats/boxplot.html
%   for your parameter estimates. Thus, we assess the quality of the MLE 
%   (via their algorithm).
 
% Crucially, you SUBTRACT THE TRUE VALUE of the parameter from your 
%   estimates before you boxplot them, so the boxplot shows the DEVIATION 
%   FROM THE TRUTH. Got it?

%clear

% input parameters
initial_df = 1; iter = 500; dim = 3;
wgts_low = 1/n_samp_low * ones(n_samp_low, 1); wgts_high = 1/n_samp_high * ones(n_samp_high, 1);
%%
% storage
% % variance-covariance matrix
%vcov_mat = zeros(dim, dim, reps);
% % random sample
%x_mat_low_store = zeros(dim, n_samp_low, reps); x_mat_high_store = zeros(dim, n_samp_high, reps);
% % nu
final_nu_low_store = zeros(reps, 1); final_nu_high_store = zeros(reps, 1);
% % mu
mu_vec_low_store = zeros(dim, reps); mu_vec_high_store = zeros(dim, reps);
% % sigma
sigma_mat_low_store = zeros(dim, dim, reps); sigma_mat_high_store = zeros(dim, dim, reps);

tic;
for r = 1:reps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get random sample of a Student t dist
    %rng(4, 'twister');
    %vcov_mat(:,:,r) = corrmat(dim);
    %x_mat_low_store(:,:,r) = mvtrnd(vcov_mat(:,:,r), true_df, n_samp_low)';
    %x_mat_high_store(:,:,r) = mvtrnd(vcov_mat(:,:,r), true_df, n_samp_high)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % call function
    [final_nu_low, nu_vec_low, mu_low, sigma_low, stop_low] = ex1a_function_MMFAlgorithm_ver3(x_mat_low_store(:,:,r), initial_df, wgts_low, iter, 0);
    [final_nu_high, nu_vec_high, mu_high, sigma_high, stop_high] = ex1a_function_MMFAlgorithm_ver3(x_mat_high_store(:,:,r), initial_df, wgts_high, iter, 0);

    % store results
    disp(['low: ', num2str(round(abs(final_nu_low - true_df),4)), ', high: ', num2str(round(abs(final_nu_high - true_df),4))]);
    if mod(r, 100) == 0; disp(['*******', num2str(r), ' reps done *******']); end
    final_nu_low_store(r) = final_nu_low; final_nu_high_store(r) = final_nu_high;
    mu_vec_low_store(:,r) = mu_low; mu_vec_high_store(:,r) = mu_high;
    sigma_mat_low_store(:,:,r) = sigma_low; sigma_mat_high_store(:,:,r) = sigma_high;
end

% report time to run
time_MMF = toc;
disp(time_MMF);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% storing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%struct_x_mat = struct('x_mat_low', x_mat_low_store, 'x_mat_high', x_mat_high_store);
%struct_final_nu = struct('final_nu_low', final_nu_low_store, 'final_nu_high', final_nu_high_store);
%struct_mu_vec = struct('mu_vec_low', mu_vec_low_store, 'mu_vec_high', mu_vec_high_store);
%struct_sigma_mat = struct('sigma_mat_low', sigma_mat_low_store, 'sigma_mat_high', sigma_mat_high_store);
%struct_comb = struct('struct_final_nu', struct_final_nu, 'struct_mu_vec', struct_mu_vec, 'struct_sigma_mat', struct_sigma_mat);
%save('ex1c_x_mat.mat', 'struct_x_mat')
%save('ex1c_params.mat', 'struct_comb')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get deviation from true values
final_nu_low_dev = final_nu_low_store - true_df; final_nu_high_dev = final_nu_high_store - true_df;
mu_vec_low_dev = mu_vec_low_store; mu_vec_high_dev = mu_vec_high_store;
sigma_mat_low_dev = sigma_mat_low_store - vcov_mat; sigma_mat_high_dev = sigma_mat_high_store - vcov_mat;

% nu
nu_plot = {final_nu_low_dev, final_nu_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(nu_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('degrees of freedom', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'df_boxplot_MMF.png')

% mu
mu_plot = {mu_vec_low_dev', mu_vec_high_dev'};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(mu_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('mean vector', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'mu_boxplot_MMF.png')

% sigma
el_sigma_tril = (size(vcov_mat(:,:,1),1) * (size(vcov_mat(:,:,1),1) + 1)) / 2;
sigma_plot_low_dev = zeros(reps, el_sigma_tril); sigma_plot_high_dev = zeros(reps, el_sigma_tril);
for r = 1:reps
    sigma_plot_low_dev(r,1) = sigma_mat_low_dev(1,1,r);
    sigma_plot_low_dev(r,2) = sigma_mat_low_dev(2,1,r);
    sigma_plot_low_dev(r,3) = sigma_mat_low_dev(2,2,r);
    sigma_plot_low_dev(r,4) = sigma_mat_low_dev(3,1,r);
    sigma_plot_low_dev(r,5) = sigma_mat_low_dev(3,2,r);
    sigma_plot_low_dev(r,6) = sigma_mat_low_dev(3,3,r);
    
    sigma_plot_high_dev(r,1) = sigma_mat_high_dev(1,1,r);
    sigma_plot_high_dev(r,2) = sigma_mat_high_dev(2,1,r);
    sigma_plot_high_dev(r,3) = sigma_mat_high_dev(2,2,r);
    sigma_plot_high_dev(r,4) = sigma_mat_high_dev(3,1,r);
    sigma_plot_high_dev(r,5) = sigma_mat_high_dev(3,2,r);
    sigma_plot_high_dev(r,6) = sigma_mat_high_dev(3,3,r);
end

sigma_plot = {sigma_plot_low_dev, sigma_plot_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(sigma_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('variance-covariance matrix', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'sigma_boxplot_MMF.png')
% 1d
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

% use same 500 replications as above as the data!!
% use x vector created above:
%    - x_mat_low_store for T = 200
%    - x_mat_high_store for T = 2000
%%
% % nu
final_nu_ll_low_store = zeros(reps, 1); final_nu_ll_high_store = zeros(reps, 1);
% % mu
mu_vec_ll_low_store = zeros(dim, reps); mu_vec_ll_high_store = zeros(dim, reps);
% % sigma
sigma_mat_ll_low_store = zeros(dim, dim, reps); sigma_mat_ll_high_store = zeros(dim, dim, reps);

tic;
for r = 1:reps
    params_low = MVTestimation(x_mat_low_store(:,:,r)');
    params_high = MVTestimation(x_mat_high_store(:,:,r)');
    sigma_low(1,1) = params_low(5); sigma_low(2,1) = params_low(6); sigma_low(3,1) = params_low(7);
    sigma_low(2,2) = params_low(8); sigma_low(3,2) = params_low(9); sigma_low(3,3) = params_low(10);
    sigma_high(1,1) = params_high(5); sigma_high(2,1) = params_high(6); sigma_high(3,1) = params_high(7);
    sigma_high(2,2) = params_high(8); sigma_high(3,2) = params_high(9); sigma_high(3,3) = params_high(10);

    disp(['low: ', num2str(round(abs(params_low(1) - true_df),4)), ', high: ', num2str(round(abs(params_high(1) - true_df),4))]);
    if mod(r, 100) == 0; disp(['*******', num2str(r), ' reps done *******']); end
    final_nu_ll_low_store(r) = params_low(1); final_nu_ll_high_store(r) = params_high(1);
    mu_vec_ll_low_store(:,r) = params_low(2:4); mu_vec_ll_high_store(:,r) = params_high(2:4);
    sigma_mat_ll_low_store(:,:,r) = sigma_low; sigma_mat_ll_high_store(:,:,r) = sigma_high;
end
time_ll = toc;
disp(time_ll);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get deviation from true values
final_nu_ll_low_dev = final_nu_ll_low_store - true_df; final_nu_ll_high_dev = final_nu_high_store - true_df;
mu_vec_ll_low_dev = mu_vec_ll_low_store; mu_vec_ll_high_dev = mu_vec_ll_high_store;
sigma_mat_ll_low_dev = sigma_mat_ll_low_store - vcov_mat; sigma_mat_ll_high_dev = sigma_mat_ll_high_store - vcov_mat;

% nu
nu_plot = {final_nu_ll_low_dev, final_nu_ll_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(nu_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('degrees of freedom', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'df_boxplot_ll.png')
% mu
mu_plot = {mu_vec_ll_low_dev', mu_vec_ll_high_dev'};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(mu_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('mean vector', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'mu_boxplot_ll.png')

% sigma
el_sigma_tril = (size(vcov_mat(:,:,1),1) * (size(vcov_mat(:,:,1),1) + 1)) / 2;
sigma_plot_ll_low_dev = zeros(reps, el_sigma_tril); sigma_plot_ll_high_dev = zeros(reps, el_sigma_tril);
for r = 1:reps
    sigma_plot_ll_low_dev(r,1) = sigma_mat_ll_low_dev(1,1,r);
    sigma_plot_ll_low_dev(r,2) = sigma_mat_ll_low_dev(2,1,r);
    sigma_plot_ll_low_dev(r,3) = sigma_mat_ll_low_dev(2,2,r);
    sigma_plot_ll_low_dev(r,4) = sigma_mat_ll_low_dev(3,1,r);
    sigma_plot_ll_low_dev(r,5) = sigma_mat_ll_low_dev(3,2,r);
    sigma_plot_ll_low_dev(r,6) = sigma_mat_ll_low_dev(3,3,r);
    
    sigma_plot_ll_high_dev(r,1) = sigma_mat_ll_high_dev(1,1,r);
    sigma_plot_ll_high_dev(r,2) = sigma_mat_ll_high_dev(2,1,r);
    sigma_plot_ll_high_dev(r,3) = sigma_mat_ll_high_dev(2,2,r);
    sigma_plot_ll_high_dev(r,4) = sigma_mat_ll_high_dev(3,1,r);
    sigma_plot_ll_high_dev(r,5) = sigma_mat_ll_high_dev(3,2,r);
    sigma_plot_ll_high_dev(r,6) = sigma_mat_ll_high_dev(3,3,r);
end

sigma_plot_ll = {sigma_plot_ll_low_dev, sigma_plot_ll_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(sigma_plot_ll, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('variance-covariance matrix', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'sigma_boxplot_ll.png')

% 1e
% Using the SAME 500 replications as above (so results are even more 
%   comparable), do the same but using the ECME algorithm from
%   C Liu and D B Rubin, (1995) "ML estimation of the t distribution using 
%   EM and its extensions, ECM and ECME", Statistica Sinica, [5, pp19-39]
%   http://www3.stat.sinica.edu.tw/statistica/oldpdf/A5n12.pdf
%   function fitt
%
%   and also the fast *approximate* method from:
%   C Aeschlimna, J Park and KA Cak, "A Novel Parameter Estimation Algorithm 
%   for the Multivariate t-Distribution and its Application to Computer 
%   Vision" [ECCV 2010]
%   http://link.springer.com/chapter/10.1007%2F978-3-642-15552-9_43
%   function fitt_approx

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st algo: ECME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% storage
% % nu
nu_ECME_low_store = zeros(reps, 1); nu_ECME_high_store = zeros(reps, 1);
% % mu
mu_vec_ECME_low_store = zeros(dim, reps); mu_vec_ECME_high_store = zeros(dim, reps);
% % sigma
sigma_mat_ECME_low_store = zeros(dim, dim, reps); sigma_mat_ECME_high_store = zeros(dim, dim, reps);

tic;
for r = 1:reps
    % call function
    [mu_ECME_low, sigma_ECME_low, nu_ECME_low] = fitt(x_mat_low_store(:,:,r)');
    [mu_ECME_high, sigma_ECME_high, nu_ECME_high] = fitt(x_mat_high_store(:,:,r)');

    % store results
    disp(['low: ', num2str(round(abs(nu_ECME_low - true_df),4)), ', high: ', num2str(round(abs(nu_ECME_high - true_df),4))]);
    if mod(r, 100) == 0; disp(['*******', num2str(r), ' reps done *******']); end
    nu_ECME_low_store(r) = nu_ECME_low; nu_ECME_high_store(r) = nu_ECME_high;
    mu_vec_ECME_low_store(:,r) = mu_ECME_low; mu_vec_ECME_high_store(:,r) = mu_ECME_high;
    sigma_mat_ECME_low_store(:,:,r) = sigma_ECME_low; sigma_mat_ECME_high_store(:,:,r) = sigma_ECME_high;
end
time_ECME = toc;
disp(time_ECME);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get deviation from true values
nu_ECME_low_dev = nu_ECME_low_store - true_df; nu_ECME_high_dev = nu_ECME_high_store - true_df;
mu_vec_ECME_low_dev = mu_vec_ECME_low_store; mu_vec_ECME_high_dev = mu_vec_ECME_high_store;
sigma_mat_ECME_low_dev = sigma_mat_ECME_low_store - vcov_mat; sigma_mat_ECME_high_dev = sigma_mat_ECME_high_store - vcov_mat;

% nu
nu_ECME_low_dev_adj = nu_ECME_low_dev .* (nu_ECME_low_dev <= 100);
nu_ECME_high_dev_adj = nu_ECME_high_dev .* (nu_ECME_high_dev <= 100);
nu_ECME_plot = {nu_ECME_low_dev_adj, nu_ECME_high_dev_adj};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(nu_ECME_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset, 'OutlierSize', 5)
title('degrees of freedom', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'df_boxplot_ECME.png')

% mu
mu_ECME_plot = {mu_vec_ECME_low_dev', mu_vec_ECME_high_dev'};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(mu_ECME_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('mean vector', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'mu_boxplot_ECME.png')

% sigma
el_sigma_tril = (size(vcov_mat(:,:,1),1) * (size(vcov_mat(:,:,1),1) + 1)) / 2;
sigma_plot_low_dev = zeros(reps, el_sigma_tril); sigma_plot_high_dev = zeros(reps, el_sigma_tril);
for r = 1:reps
    sigma_plot_low_dev(r,1) = sigma_mat_ECME_low_dev(1,1,r);
    sigma_plot_low_dev(r,2) = sigma_mat_ECME_low_dev(2,1,r);
    sigma_plot_low_dev(r,3) = sigma_mat_ECME_low_dev(2,2,r);
    sigma_plot_low_dev(r,4) = sigma_mat_ECME_low_dev(3,1,r);
    sigma_plot_low_dev(r,5) = sigma_mat_ECME_low_dev(3,2,r);
    sigma_plot_low_dev(r,6) = sigma_mat_ECME_low_dev(3,3,r);
    
    sigma_plot_high_dev(r,1) = sigma_mat_ECME_high_dev(1,1,r);
    sigma_plot_high_dev(r,2) = sigma_mat_ECME_high_dev(2,1,r);
    sigma_plot_high_dev(r,3) = sigma_mat_ECME_high_dev(2,2,r);
    sigma_plot_high_dev(r,4) = sigma_mat_ECME_high_dev(3,1,r);
    sigma_plot_high_dev(r,5) = sigma_mat_ECME_high_dev(3,2,r);
    sigma_plot_high_dev(r,6) = sigma_mat_ECME_high_dev(3,3,r);
end

sigma_plot = {sigma_plot_low_dev, sigma_plot_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(sigma_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('variance-covariance matrix', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'sigma_boxplot_ECME.png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd algo: approx %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note that c is the median not the mean of the data
%%
% storage
% % nu
nu_approx_low_store = zeros(reps, 1); nu_approx_high_store = zeros(reps, 1);
% % mu
mu_vec_approx_low_store = zeros(dim, reps); mu_vec_approx_high_store = zeros(dim, reps);
% % sigma
sigma_mat_approx_low_store = zeros(dim, dim, reps); sigma_mat_approx_high_store = zeros(dim, dim, reps);

tic;
for r = 1:reps
    % call function
    [mu_approx_low, sigma_approx_low, nu_approx_low] = fitt_approx(x_mat_low_store(:,:,r)');
    [mu_approx_high, sigma_approx_high, nu_approx_high] = fitt_approx(x_mat_high_store(:,:,r)');

    % store results
    disp(['low: ', num2str(round(abs(nu_approx_low - true_df),4)), ', high: ', num2str(round(abs(nu_approx_high - true_df),4))]);
    if mod(r, 100) == 0; disp(['*******', num2str(r), ' reps done *******']); end
    nu_approx_low_store(r) = nu_approx_low; nu_approx_high_store(r) = nu_approx_high;
    mu_vec_approx_low_store(:,r) = mu_approx_low; mu_vec_approx_high_store(:,r) = mu_approx_high;
    sigma_mat_approx_low_store(:,:,r) = sigma_approx_low; sigma_mat_approx_high_store(:,:,r) = sigma_approx_high;
end
time_approx = toc;
disp(time_approx);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get deviation from true values
nu_approx_low_dev = nu_approx_low_store - true_df; nu_approx_high_dev = nu_approx_high_store - true_df;
mu_vec_approx_low_dev = mu_vec_approx_low_store; mu_vec_approx_high_dev = mu_vec_approx_high_store;
sigma_mat_approx_low_dev = sigma_mat_approx_low_store - vcov_mat; sigma_mat_approx_high_dev = sigma_mat_approx_high_store - vcov_mat;

% nu
nu_approx_plot = {nu_approx_low_dev, nu_approx_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(nu_approx_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset, 'OutlierSize', 5)
title('degrees of freedom', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'df_boxplot_approx.png')

% mu
mu_approx_plot = {mu_vec_approx_low_dev', mu_vec_approx_high_dev'};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(mu_approx_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('mean vector', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'mu_boxplot_approx.png')

% sigma
el_sigma_tril = (size(vcov_mat(:,:,1),1) * (size(vcov_mat(:,:,1),1) + 1)) / 2;
sigma_plot_low_dev = zeros(reps, el_sigma_tril); sigma_plot_high_dev = zeros(reps, el_sigma_tril);
for r = 1:reps
    sigma_plot_low_dev(r,1) = sigma_mat_low_dev(1,1,r);
    sigma_plot_low_dev(r,2) = sigma_mat_low_dev(2,1,r);
    sigma_plot_low_dev(r,3) = sigma_mat_low_dev(2,2,r);
    sigma_plot_low_dev(r,4) = sigma_mat_low_dev(3,1,r);
    sigma_plot_low_dev(r,5) = sigma_mat_low_dev(3,2,r);
    sigma_plot_low_dev(r,6) = sigma_mat_low_dev(3,3,r);
    
    sigma_plot_high_dev(r,1) = sigma_mat_high_dev(1,1,r);
    sigma_plot_high_dev(r,2) = sigma_mat_high_dev(2,1,r);
    sigma_plot_high_dev(r,3) = sigma_mat_high_dev(2,2,r);
    sigma_plot_high_dev(r,4) = sigma_mat_high_dev(3,1,r);
    sigma_plot_high_dev(r,5) = sigma_mat_high_dev(3,2,r);
    sigma_plot_high_dev(r,6) = sigma_mat_high_dev(3,3,r);
end

sigma_plot = {sigma_plot_low_dev, sigma_plot_high_dev};
labelset = {'A', 'B'};

figure('Position', [400 75 500 300])
boxplotGroup(sigma_plot, 'interGroupSpace', 5, 'PrimaryLabels', labelset)
title('variance-covariance matrix', 'FontName', 'FixedWidth')
subtitle(['A: sample size = ', num2str(n_samp_low), '  B: sample size: ', num2str(n_samp_high)])
saveas(gcf,'sigma_boxplot_approx.png')
