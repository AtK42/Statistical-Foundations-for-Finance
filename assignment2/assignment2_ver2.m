%% exercise 1 (from first mail, 26.10.2022)
% define parameters
delim = '************************************';
reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
df = 2; loc = 1; scale = 2; alpha = .1;
initvec = [2 2 2]; % prior assumption for mle: [df, location, scale]

% set seed
rng(6, 'twister')

% initialize variables
% % matlab
ES_vec = zeros(n_BS, 1);
ci_length_para = zeros([reps length(n_samp_vec)]);
coverage_para = zeros([reps length(n_samp_vec)]);
% % MP
ES_vec_MP = zeros(n_BS, 1);
ci_length_para_MP = zeros([reps length(n_samp_vec)]);
coverage_para_MP = zeros([reps length(n_samp_vec)]);

% calculate theoretical ES for the student t (analytically)
c01=tinv(alpha , df);
cLS = loc+scale*c01; % cLS is cutoff Location Scale
ES_01_analytic = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1);
ES_LS_analytic = loc+scale*ES_01_analytic;
trueES = ES_LS_analytic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonparametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(delim); disp('non-parametric bootstrap');
[average_length, coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 1, [scale, loc, df], trueES, alpha);
disp(['Average nonparametric CI Length: ', num2str(average_length, '%.4f')]);
disp(['Nonparametric Coverage Ratio:    ', num2str(coverage_ratio, '%.4f')]);

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(delim); disp('parametric bootstrap');
for k = 1:length(n_samp_vec)
    disp(delim); disp(['calculations for sample size = ', num2str(k), ' in progress']);
    for i = 1:reps
        
        % generate random sample of a (regular) loc-scale t dist
        data = loc + scale * trnd(df, n_samp_vec(k), 1);
       
        for j = 1:n_BS
            
            % create bootstrap sample
            ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
            bs_samp = data(ind);
            
            % parametric
            % % via matlab's mle function
            para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
            para_loc_hat = para_bs_hat(1); para_scale_hat = para_bs_hat(2); para_df_hat = para_bs_hat(3);
            
            % % via MP function from his book
            para_bs_hat_MP = tlikmax(bs_samp, initvec);
            para_loc_hat_MP = para_bs_hat_MP(1); para_scale_hat_MP = para_bs_hat_MP(2); para_df_hat_MP = para_bs_hat_MP(3);

            % calculate theoretical ES based on the parameter estimates
            % % matlab
            c01 = tinv(alpha , para_df_hat);
            ES_vec(j) = para_loc_hat + para_scale_hat * (-tpdf(c01,para_df_hat)/tcdf(c01,para_df_hat) * (para_df_hat+c01^2)/(para_df_hat-1));
            % % MP
            c01_MP = tinv(alpha , para_df_hat_MP);
            ES_vec_MP(j) = para_loc_hat_MP + para_scale_hat_MP * (-tpdf(c01,para_df_hat_MP)/tcdf(c01,para_df_hat_MP) * (para_df_hat_MP+c01^2)/(para_df_hat_MP-1));
        end % j-loop (bootstrap loop)

        % compute length of the CI and coverage
        % % matlab
        ci_para = quantile(ES_vec, [alpha/2 1-alpha/2]);
        low_para = ci_para(1); high_para = ci_para(2);
        ci_length_para(i, k) = high_para - low_para;
        if ES_LS_analytic >= low_para && ES_LS_analytic <= high_para
            coverage_para(i, k) = 1;
        end
        % % MP
        ci_para_MP = quantile(ES_vec_MP, [alpha/2 1-alpha/2]);
        low_para_MP = ci_para_MP(1); high_para_MP = ci_para_MP(2);
        ci_length_para_MP(i, k) = high_para_MP - low_para_MP;
        if ES_LS_analytic >= low_para_MP && ES_LS_analytic <= high_para_MP
            coverage_para_MP(i, k) = 1;
        end
        
        if mod(i, 10) == 0
           disp(['rep ', num2str(i), ' out of ', reps, ' done']);
        end
    end %i-loop
end %k-loop

% compare the length of the CIs as a function of the sample size n_samp (= T)
disp(delim); disp('CI length using'); disp(delim)
disp(["matlab function: ", num2str(mean(ci_length_para),    '%.4f')]);
disp(["MP's function:   ", num2str(mean(ci_length_para_MP), '%.4f')]);

% compare coverage ratio of the theoretical ES as a function of the sample size n_samp (= T)
disp(delim); disp('coverage ratio using'); disp(delim)
disp(["matlab function: ", num2str(mean(coverage_para),    '%.4f')]);
disp(["MP's function:   ", num2str(mean(coverage_para_MP), '%.4f')]);

disp(delim);
%% exercise 2
% define parameters
delim = '************************************';
n_samp = 1e7; n_samp_vec = [250 500 2000];
df_vec = [3 6]; % degrees of freedom of the NCT
mu_vec = [-3 -2 -1 0]; % (numerator) non-centrality parameter of the NCT
theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
seed = 6; alpha = .1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 1 (learn how to simulate from non-central location-scale t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see function asymtrnd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 2 (calculate true ES of the NTC for the following parameters and via
%   (i) simulation and 
%   (ii) integral definition of the NTC for):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different df:                           |    3|    6|     |     |
% % four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|

disp(delim); disp(['ES for different df and non-centrality', newline, 'parameters of the NCT (mu)']);
for df = 1:numel(df_vec)
        disp(delim); disp(['for df = ', num2str(df_vec(df))]);
    % (i) Simulation:
        ES_sim = zeros(numel(mu_vec), 1);
        for mu = 1:numel(mu_vec)
                ES_sim(mu) = Simulated_ES_NCT(n_samp, df_vec(df), mu_vec(mu), alpha, seed);
        end % end mu-loop
        disp(['via Simulation:          ', num2str(ES_sim', ' %.4f')]);

    % (ii) numeric integration:
        ES_num = zeros(numel(mu_vec), 1);
        for mu = 1:numel(mu_vec)
                c01=nctinv(alpha , df_vec(df), mu_vec(mu));
                I01 = @(x) x.*nctpdf(x, df_vec(df), mu_vec(mu)); %note that the problem with nctpdf mentioned in footnote 11 on p.373 in the intermediate prob book has been solved in the standard matlab function, hence it is used here
                ES_num(mu) = integral(I01 , -Inf , c01) / alpha;
        end % end mu-loop
        disp(['via Numeric Integration: ', num2str(ES_num', '  %.4f')]);
end % df-loop (for df)
disp(delim);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 3 (report average length of CI and actual coverage with both bootstrap methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different df:                           |    3|    6|     |     |
% % three different sample sizes (n_samp):      |  250|  500| 2000|     |
% % four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|
delim = '************************************';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for df = 1:numel(df_vec)
    disp(delim); disp(['for df = ', df]);
    for mu=1:numel(mu_vec)
        [average_length, coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 2, [df_vec(df), mu_vec(mu)], ES_num(mu), alpha);
        disp(['Average nonparametric CI Length with mu = ', num2str(mu_vec(mu)),  ': ', num2str(average_length, '%.4f')]);
        disp(['Nonparametric Coverage Ratio: with mu =   ', num2str(mu_vec(mu)),  ': ', num2str(coverage_ratio, '%.4f')]);
    end % mu-loop
end

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
%!!Coverage rate is wrong!!
for df = 1:numel(df_vec)
    disp(delim); disp(['for df = ', df]);
    for k = 1:length(n_samp_vec)
        for mu = 1:numel(mu_vec)
            for i = 1:reps
                % generate random sample of a (regular) loc-scale t dist
                data = loc + scale * asymtrnd(n_samp_vec(k), mu_vec(mu), df_vec(df), 1);
                
                for j = 1:n_BS
                    bs_samp = para_loc_hat + para_scale_hat * trnd(para_df_hat, n_samp_vec(k), 1);
                    VaR = quantile(bs_samp, alpha);
                    temp = bs_samp(bs_samp<=VaR);
                    ES_vec(j) = mean(temp);
                end
                ci_para = quantile(ES_vec, [alpha/2 1-alpha/2]);
                low_para = ci_para(1); high_para = ci_para(2);
                ci_length_para(i, mu) = high_para - low_para;
                if ES_LS_analytic >= low_para && ES_LS_analytic <= high_para
                    coverage_para(i,mu) = 1;
                end
            end % end i
        end % end mu
        disp(mean(ci_length_para));
        disp(mean(coverage_para));
    end % end k
end % end df

%% exercise 3









