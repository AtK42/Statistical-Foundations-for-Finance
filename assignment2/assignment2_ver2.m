%% exercise 1 (from first mail, 26.10.2022)
tic
% define parameters
delim = '************************************';
reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
df = 4; loc = 1; scale = 2; alpha = .01; a = 0.1;
para_method = "MP";

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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonparametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(delim); disp('non-parametric bootstrap');
[ci_length_nonpara, coverage_ratio_nonpara, average_length_nonpara, mean_coverage_ratio_nonpara] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 1, [scale, loc, df], trueES, alpha, a);
disp('***');
disp(['Average nonparametric CI Length: ', num2str(average_length_nonpara, '% 7.4f')]);
disp(['Nonparametric Coverage Ratio:    ', num2str(mean_coverage_ratio_nonpara, '% 7.4f')]);

struct_nonpara = struct('average_length', average_length_nonpara, 'ci_length', ci_length_nonpara, 'mean_coverage_ratio', mean_coverage_ratio_nonpara, 'coverage_ratio', coverage_ratio_nonpara);

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
%ci_length = zeros([reps 1]);
%coverage = zeros([reps 1]);

disp(delim); disp(delim); disp('parametric bootstrap');
for k = 1:length(n_samp_vec)
    disp(delim); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for i = 1:reps
        
        % generate random sample of a (regular) loc-scale t dist
        data = loc + scale * trnd(df, n_samp_vec(k), 1);
        initvec = [df loc scale]; % prior assumption for mle: [df, location, scale]
        para_bs_hat = tlikmax(data, initvec);
        para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

        sim_para_mat = para_loc_hat + para_scale_hat * trnd(para_df_hat, [n_samp_vec(k) n_BS]);
        ES_vec = zeros([1 n_BS]);
        for j = 1:n_BS
            % create bootstrap sample
            temp = sim_para_mat(:, j);
            params_temp = mle(temp, 'Distribution', 'tLocationScale'); %output [loc scale df]
            params_temp_MP = tlikmax(temp, initvec); %[df loc scale]
            %ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
            %bs_samp = data(ind);
            
            % parametric
            % % via matlab's mle function
            %para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
            %para_loc_hat = para_bs_hat(1); para_scale_hat = para_bs_hat(2); para_df_hat = para_bs_hat(3);
            
            % % via MP function from his book
            %para_bs_hat_MP = tlikmax(bs_samp, initvec);
            %para_df_hat_MP = para_bs_hat_MP(1); para_loc_hat_MP = para_bs_hat_MP(2); para_scale_hat_MP = para_bs_hat_MP(3);

            % calculate theoretical ES based on the parameter estimates
            % % matlab
            c01 = tinv(alpha , params_temp(3));
            ES_vec(j) = params_temp(1) + params_temp(2) * (-tpdf(c01,params_temp(3))/tcdf(c01,params_temp(3)) * (params_temp(3)+c01^2)/(params_temp(3)-1));
            % % MP
            c01_MP = tinv(alpha , params_temp_MP(1));
            ES_vec_MP(j) = params_temp_MP(2) + params_temp_MP(3) * (-tpdf(c01,params_temp_MP(1))/tcdf(c01,params_temp_MP(1)) * (params_temp_MP(1)+c01^2)/(params_temp_MP(1)-1));
        end % j-loop (bootstrap loop)

        % compute length of the CI and coverage
        % % matlab
        ci_para = quantile(ES_vec, [a/2 1-a/2]);
        disp(ci_para);
        low_para = ci_para(1); high_para = ci_para(2);
        ci_length_para(i, k) = high_para - low_para;
        if ES_LS_analytic >= low_para && ES_LS_analytic <= high_para
            coverage_para(i, k) = 1;
        end
        
        % % MP
        ci_para_MP = quantile(ES_vec_MP, [a/2 1-a/2]);
        low_para_MP = ci_para_MP(1); high_para_MP = ci_para_MP(2);
        ci_length_para_MP(i, k) = high_para_MP - low_para_MP;
        if ES_LS_analytic >= low_para_MP && ES_LS_analytic <= high_para_MP
            coverage_para_MP(i, k) = 1;
        end
        
        if mod(i, 10) == 0
           disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
        end
    end %i-loop
end %k-loop

disp('***');
% compare the length of the CIs as a function of the sample size n_samp (= T)
disp(delim); disp('CI length using');
disp(["matlab function: ", num2str(mean(ci_length_para),    '% 7.4f')]);
disp(["MP's function:   ", num2str(mean(ci_length_para_MP), '% 7.4f')]);

% compare coverage ratio of the theoretical ES as a function of the sample size n_samp (= T)
disp(delim); disp('coverage ratio using');
disp(["matlab function: ", num2str(mean(coverage_para),    '% 7.4f')]);
disp(["MP's function:   ", num2str(mean(coverage_para_MP), '% 7.4f')]);

disp(delim);

% save output
struct_para = struct('mean_ci_length_para'   , mean(ci_length_para)   , 'ci_length_para'   , ci_length_para   , 'mean_coverage_para'   , mean(coverage_para)   , 'coverage_para'   , coverage_para);
struct_para_MP = struct('mean_ci_length_para_MP', mean(ci_length_para_MP), 'ci_length_para_MP', ci_length_para_MP, 'mean_coverage_para_MP', mean(coverage_para_MP), 'coverage_para_MP', coverage_para_MP);
struct_comb = struct('struct_nonpara', struct_nonpara, 'struct_para', struct_para, 'struct_para_MP', struct_para_MP);
save('ex1.mat', 'struct_comb')

% struct
%         | n_samp = 250 | n_samp = 500 | n_samp = 2000 |
%   ...   |              |              |               |

toc
disp(delim); disp(delim); %took 16267.616772 seconds
%% exercise 2 (simulate from NCT (part 1) & calculate true ES of NCT (part 2))
tic
% define parameters
delim = '************************************';
n_samp = 1e7; loc = 1; scale = 2;
%reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
df_vec = [3 6]; % degrees of freedom of the NCT
mu_vec = [-2 -1 0]; % (numerator) non-centrality parameter of the NCT
%theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
seed = rand*1000; alpha = .01; a = 0.1;

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

ES_sim = zeros(numel(mu_vec), numel(df_vec));
ES_num = zeros(numel(mu_vec), numel(df_vec));

disp(delim); disp(['ES for different df and non-centrality', newline, 'parameters of the NCT (mu)']);
for df = 1:numel(df_vec)
        disp(delim); disp(['for df = ', num2str(df_vec(df))]);
    % (i) Simulation:
        for mu = 1:numel(mu_vec)
                ES_sim(mu, df) = loc + scale * Simulated_ES_NCT(n_samp, df_vec(df), mu_vec(mu), alpha, seed);
        end % end mu-loop
        disp(['via Simulation:          ', num2str(ES_sim(:, df)', '% 7.4f')]);

    % (ii) numeric integration:
        for mu = 1:numel(mu_vec)
                c01=nctinv(alpha , df_vec(df), mu_vec(mu));
                I01 = @(x) x.*nctpdf(x, df_vec(df), mu_vec(mu)); %note that the problem with nctpdf mentioned in footnote 11 on p.373 in the intermediate prob book has been solved in the standard matlab function, hence it is used here
                ES_num(mu, df) = loc + scale * integral(I01 , -Inf , c01) / alpha;
        end % end mu-loop
        disp(['via Numeric Integration: ', num2str(ES_num(:, df)', '% 7.4f')]);
end % df-loop (for df)

disp(delim);
struct_ES = struct('ES_sim', ES_sim, 'ES_num', ES_num);
save('ex2_trueES.mat', 'struct_ES');
% struct result:
%         | df = 3 | df = 6 |
% mu = -3 |        |        |
% mu = -2 |        |        |
% mu = -1 |        |        |
% mu = 0  |        |        |
toc
disp(delim); disp(delim);  %took 20.562306 seconds
%% exercise 2 (bootstrapping for first df (part 3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 3 - first df (report average length of CI and actual coverage with both bootstrap methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different df:                           |    3|    6|     |     |
% % three different sample sizes (n_samp):      |  250|  500| 2000|     |
% % four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|
tic
delim = '************************************';
loc = 1; scale = 2;
reps = 200; n_samp_vec = [250, 500, 2000]; n_BS = 1000; % note that n_samp = T
df = 3; % degrees of freedom of the NCT
n_df = 1;
mu_vec = [-3 -2 -1 0]; % (numerator) non-centrality parameter of the NCT
%theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
% seed = rand*1000;
alpha = .01; a = 0.1;

ci_length = zeros([reps 1]);
coverage = zeros([reps 1]);

%ci_length_nonpara = zeros(numel(mu), reps, numel(n_samp_vec));
%coverage_ratio_nonpara = zeros(numel(mu), reps, numel(n_samp_vec));
%ci_length_para = zeros(numel(mu), reps, numel(n_samp_vec));
%coverage_para = zeros(numel(mu), reps, numel(n_samp_vec));

ci_length_nonpara = cat(4, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4
coverage_ratio_nonpara = cat(4, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4
ci_length_para = cat(4, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4
coverage_para = cat(4, ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(['for df = ', num2str(df)]);
for mu=1:numel(mu_vec)
    [ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 2, [scale, loc, df, mu_vec(mu)], ES_num(mu, n_df), alpha, a);
    disp(['Average nonparametric CI Length with mu = ', num2str(mu_vec(mu)),  ': ', num2str(average_length, '% 7.4f')]);
    disp(['Nonparametric Coverage Ratio: with mu =   ', num2str(mu_vec(mu)),  ': ', num2str(mean_coverage_ratio, '% 7.4f')]);
    ci_length_nonpara(:, :, mu) = ci_length;
    coverage_ratio_nonpara(:, :, mu) = coverage_ratio;
end % mu-loop

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
%!!check coverage rate (might be wrong still)!!
disp(delim); disp(['ex2: parametric bootstrap for df = ', num2str(df)]);
for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for mu = 1:numel(mu_vec)
        disp(['   starting calculations for non-centrality param mu = ', num2str(mu_vec(mu))]);
        ci_length = zeros([reps 1]);
        coverage = zeros([reps 1]);

        for i = 1:reps

            % generate random sample of a (regular) loc-scale t dist
            data = loc + scale * asymtrnd([n_samp_vec(k) 1], mu_vec(mu), df);
            initvec = [df loc scale]; % [df loc scale]
            para_bs_hat = tlikmax(data, initvec);
            para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);
            
            sim_para_mat = para_loc_hat + para_scale_hat * trnd(para_df_hat, [n_samp_vec(k) n_BS]);
            ES_vec = zeros([1 n_BS]);
            for j = 1:n_BS
                temp = sim_para_mat(:, j);
                initvec = [para_df_hat, para_loc_hat, para_scale_hat];
                params_temp = tlikmax(temp, initvec); %[df, loc, scale]

                % create bootstrap sample
                %ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
                %bs_samp = data(ind);

                % parametric
                %para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
                %para_loc_hat = para_bs_hat(1); para_scale_hat = para_bs_hat(2); para_df_hat = para_bs_hat(3);
                %para_bs_hat = tlikmax(bs_samp, initvec);
                %para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

                % calculate theoretical ES based on the parameter estimates
                %c01 = tinv(alpha , para_df_hat);
                %ES_vec(j) = para_loc_hat + para_scale_hat * (-tpdf(c01,para_df_hat)/tcdf(c01,para_df_hat) * (para_df_hat+c01^2)/(para_df_hat-1));

                c01 = tinv(alpha, params_temp(1));
                ES_vec(j) = params_temp(2) + params_temp(3) * (-tpdf(c01, params_temp(1))/tcdf(c01, params_temp(1)) * (params_temp(1)+c01^2)/(params_temp(1)-1));
            end % j-loop

            % compute length of the CI and coverage
            ci_para = quantile(ES_vec, [a/2 1-a/2]);
            low_para = ci_para(1); high_para = ci_para(2);
            %ci_length_para(mu, i, k) = high_para - low_para;
            ci_length(i) = high_para - low_para;
            if ES_num(mu, n_df) >= low_para && ES_num(mu, n_df) <= high_para
                %coverage_para(mu, i, k) = 1;
                coverage(i) = 1;
            end

            if mod(i, 10) == 0
                disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
            end
        end % i-loop (reps)
        ci_length_para(:, k, 1, mu) = ci_length;
        coverage_para(:, k, 1, mu) = coverage;
        % ci_length_para(:, k, 1, mu) = ci_length(:, k);
        % coverage_para(:, k, 1, mu) = coverage(:, k);
    end % mu-loop

    %disp(mean(ci_length_para));
    %disp(mean(coverage_para));
end % k-loop (samp size)

% save
struct_nonpara_firstdf = struct('average_length', mean(ci_length_nonpara), 'ci_length', ci_length_nonpara, 'mean_coverage_ratio', mean(coverage_ratio_nonpara), 'coverage_ratio', coverage_ratio_nonpara);
struct_para_firstdf = struct('mean_ci_length_para', mean(ci_length_para), 'ci_length_para', ci_length_para, 'mean_coverage_ratio_para', mean(coverage_para), 'coverage_ratio_para', coverage_para);
struct_comb = struct('struct_nonpara_firstdf', struct_nonpara_firstdf, 'struct_para_firstdf', struct_para_firstdf);
% save('results/ex2_firstdf_len+coverage.mat', 'struct_comb');
save('ex2_firstdf_len+coverage.mat', 'struct_comb');

toc
disp(delim); disp(delim);  %took 22187.957696 seconds
%% exercise 2 (bootstrapping for second df (part 3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 3 - second df (report average length of CI and actual coverage with both bootstrap methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different df:                           |    3|    6|     |     |
% % three different sample sizes (n_samp):      |  250|  500| 2000|     |
% % four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|
tic
delim = '************************************';
loc = 1; scale = 2;
reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
df = 6; % degrees of freedom of the NCT
n_df = 2;
mu_vec = [-3 -2 -1 0]; % (numerator) non-centrality parameter of the NCT
%theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
%seed = rand*1000;
alpha = .01; a = 0.1;

ci_length = zeros([reps numel(n_samp_vec)]);
coverage = zeros([reps numel(n_samp_vec)]);

%ci_length_nonpara = zeros(numel(mu), reps, numel(n_samp_vec));
%coverage_ratio_nonpara = zeros(numel(mu), reps, numel(n_samp_vec));
%ci_length_para = zeros(numel(mu), reps, numel(n_samp_vec));
%coverage_para = zeros(numel(mu), reps, numel(n_samp_vec));

ci_length_nonpara = cat(4, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4
coverage_ratio_nonpara = cat(4, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4
ci_length_para = cat(4, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4
coverage_para = cat(4, ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(['ex2: non-parametric bootstrap for df = ', num2str(df)]);
for mu=1:numel(mu_vec)
    [ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 2, [scale, loc, df, mu_vec(mu)], ES_num(mu, n_df), alpha, a);
    disp(['Average nonparametric CI Length with mu = ', num2str(mu_vec(mu)),  ': ', num2str(average_length, '% 7.4f')]);
    disp(['Nonparametric Coverage Ratio: with mu =   ', num2str(mu_vec(mu)),  ': ', num2str(mean_coverage_ratio, '% 7.4f')]);
    ci_length_nonpara(:, :, mu) = ci_length;
    coverage_ratio_nonpara(:, :, mu) = coverage_ratio;
end % mu-loop

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%

%!!check coverage rate (might be wrong still)!!
disp(delim); disp(['ex2: parametric bootstrap for df = ', num2str(df)]);
for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for mu = 1:numel(mu_vec)
    disp(['   starting calculations for non-centrality param mu = ', num2str(mu_vec(mu))]);
        ci_length = zeros([reps 1]);
        coverage = zeros([reps 1]);

        for i = 1:reps
            % generate random sample of a (regular) loc-scale t dist
            data = loc + scale * asymtrnd([n_samp_vec(k) 1], mu_vec(mu), df);
            initvec = [df loc scale]; % [df loc scale]
            para_bs_hat = tlikmax(data, initvec);
            para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

            sim_para_mat = para_loc_hat + para_scale_hat * trnd(para_df_hat, [n_samp_vec(k) n_BS]);
            ES_vec = zeros([1 n_BS]);
            for j = 1:n_BS
                temp = sim_para_mat(:, j);
                params_temp = tlikmax(temp, initvec); %[df, loc, scale]

                % create bootstrap sample
                %ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
                %bs_samp = data(ind);

                % parametric
                %para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
                %para_loc_hat = para_bs_hat(1); para_scale_hat = para_bs_hat(2); para_df_hat = para_bs_hat(3);

                %para_bs_hat = tlikmax(bs_samp, initvec);
                %para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

                % calculate theoretical ES based on the parameter estimates
                c01 = tinv(alpha, params_temp(1));
                ES_vec(j) = params_temp(2) + params_temp(3) * (-tpdf(c01, params_temp(1))/tcdf(c01, params_temp(1)) * (params_temp(1)+c01^2)/(params_temp(1)-1));

            end % j-loop

            % compute length of the CI and coverage
            ci_para = quantile(ES_vec, [a/2 1-a/2]);
            low_para = ci_para(1); high_para = ci_para(2);
            %ci_length_para(mu, i, k) = high_para - low_para;
            ci_length(i) = high_para - low_para;
            if ES_num(mu, n_df) >= low_para && ES_num(mu, n_df) <= high_para
                %coverage_para(mu, i, k) = 1;
                coverage(i) = 1;
            end
            
            if mod(i, 10) == 0
                disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
            end
        end % i-loop (reps)

        ci_length_para(:, k, 1, mu) = ci_length;
        coverage_para(:, k, 1, mu) = coverage;
    end % mu-loop
    
    %disp(mean(ci_length_para));
    %disp(mean(coverage_para));
end % k-loop (samp size)

% save
struct_nonpara_seconddf = struct('average_length', mean(ci_length_nonpara), 'ci_length', ci_length_nonpara, 'mean_coverage_ratio', mean(coverage_ratio_nonpara), 'coverage_ratio', coverage_ratio_nonpara);
struct_para_seconddf = struct('mean_ci_length_para', mean(ci_length_para), 'ci_length_para', ci_length_para, 'mean_coverage_ratio_para', mean(coverage_para), 'coverage_ratio_para', coverage_para);
struct_comb = struct('struct_nonpara_seconddf', struct_nonpara_seconddf, 'struct_para_seconddf', struct_para_seconddf);

save('ex2_seconddf_len+coverage.mat', 'struct_comb');

toc
disp(delim); disp(delim);  %took 20413.9 seconds
%% exercise 3 (calculate true ES of stable)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate true ES of the stable dist for the following parameters and via
%   (i) simulation and 
%   (ii) Stoyanov (integral definition of the stable) for:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different tail_indices:                 |    1.6|    1.8|       |       |
tic
delim = '************************************';
loc = 1; scale = 2;
tail_index_vec = [1.6, 1.8];
n_samp = 1e8;
%reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
seed = rand*1000; alpha = .01; a = 0.1;

ES_sim = zeros(numel(tail_index_vec), 1); % mümmer glaub ich ned mache 
ES_stoy = zeros(numel(tail_index_vec), 1);

disp(delim); disp(['ES for different tail indices', newline, 'of the (symmetric) stable distribution']); disp('***');
for i = 1:numel(tail_index_vec)

    disp(['calculating for tail_index = ', num2str(tail_index_vec(i))]);

    % (i) Simulation:
    ES_sim(i) = Simulated_ES_symStable(n_samp, tail_index_vec(i), scale, loc, alpha, seed);

    % (ii) Stoyanov
    ES_stoy(i) = asymstableES(alpha , tail_index_vec(i), 0, loc, scale, 1);
end % i-loop (tail_index)

disp(delim);
disp(['via Simulation: ', num2str(ES_sim', ' % .4f')]);
disp(['via Stoyanov:   ', num2str(ES_stoy', ' % .4f')]);

disp(delim);
struct_ES = struct('ES_stoy', ES_stoy, 'ES_sim', ES_sim);
save('ex3_trueES.mat', 'struct_ES');
% struct result:
%                  |  ...   |
% tail_index = 1.6 |        |
% tail_index = 1.8 |        |

toc
disp(delim); disp(delim);  %took seconds
%% exercise 3 (bootstrapping for first tail_index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 2 - (report average length of CI and actual coverage with both bootstrap methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different tail_indices:                 |    1.6|    1.8|       |       |
tic
delim = '************************************';
loc = 1; scale = 2;
tail_index = 1.6; n_tail_index = 1;
reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
%seed = rand*1000;
alpha = .01; a = 0.1;
dist = 3; % corresponds to the (symmetric) stable dist

ci_length_para = zeros([reps length(n_samp_vec)]);
coverage_para = zeros([reps length(n_samp_vec)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(['ex3: non-parametric bootstrap for tail_index = ', num2str(tail_index)]);
params = [tail_index, scale, loc];
[ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, dist, params, ES_stoy(n_tail_index), alpha, a);
disp(['Average nonparametric CI Length: ', num2str(average_length, ' % 7.4f')]);
disp(['Nonparametric Coverage Ratio:    ', num2str(mean_coverage_ratio, ' % 7.4f')]);

struct_nonpara_firsttailindex = struct('average_length', average_length, 'ci_length', ci_length, 'mean_coverage_ratio', mean_coverage_ratio, 'coverage_ratio', coverage_ratio);

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
ci_length = zeros([reps 1]);
coverage = zeros([reps 1]);

disp(delim); disp(['ex3: parametric bootstrap for tail_index = ', num2str(tail_index)]);
for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for i = 1:reps
        % generate random sample from a (symmetric) stable dist
        data = stabgen(n_samp_vec(k), tail_index, 0, scale ,loc); %don't set the seed, function takes random one
        initvec = [3 loc scale];

        para_bs_hat = tlikmax(data, initvec);
        para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);
    
        sim_para_mat = para_loc_hat + para_scale_hat * trnd(para_df_hat, [n_samp_vec(k) n_BS]);
        ES_vec = zeros([1 n_BS]);
        for j = 1:n_BS
            temp = sim_para_mat(:, j);
            params_temp = tlikmax(temp, initvec); %[df, loc, scale]
            % create bootstrap sample
            %ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
            %bs_samp = data(ind);

            % parametric
            %para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
            %para_loc_hat = para_bs_hat(1); para_scale_hat = para_bs_hat(2); para_df_hat = para_bs_hat(3);

            %para_bs_hat = tlikmax(bs_samp, initvec);
            %para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

            % calculate theoretical ES based on the parameter estimates
            c01 = tinv(alpha, params_temp(1));
            ES_vec(j) = params_temp(2) + params_temp(3) * (-tpdf(c01, params_temp(1))/tcdf(c01, params_temp(1)) * (params_temp(1)+c01^2)/(params_temp(1)-1));
        end % j-loop

        % compute length of the CI and coverage
        ci_para = quantile(ES_vec, [a/2 1-a/2]);
        low_para = ci_para(1); high_para = ci_para(2);
        ci_length_para(i, k) = high_para - low_para;
        if ES_stoy(n_tail_index) >= low_para && ES_stoy(n_tail_index) <= high_para
            coverage_para(i, k) = 1;
        end
        
        if mod(i, 10) == 0
            disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
        end
    end % i-loop (reps)
    disp(k);
end % k-loop (sample size)
disp(['Average parametric CI Length: ', num2str(mean(ci_length_para), ' % 7.4f')]);
disp(['parametric Coverage Ratio:    ', num2str(mean(coverage_para), ' % 7.4f')]);

%save 
struct_para_firsttailindex = struct('mean_ci_length_para', mean(ci_length_para), 'ci_length_para', ci_length_para, 'mean_coverage_ratio_para', mean(coverage_para), 'coverage_ratio_para', coverage_para);
struct_comb = struct('struct_nonpara_firsttailindex', struct_nonpara_firsttailindex, 'struct_para_firsttailindex', struct_para_firsttailindex);
save('ex3_firsttailindex_len+coverage.mat', 'struct_comb');

toc
disp(delim); disp(delim);  %took seconds
%% exercise 3 (bootstrapping for second tail_index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 2 - (report average length of CI and actual coverage with both bootstrap methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different tail_indices:                 |    1.6|    1.8|       |       |
tic
delim = '************************************';
loc = 1; scale = 2;
tail_index = 1.8; n_tail_index = 2;
reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
%seed = rand*1000;
alpha = .01; a = 0.1;
dist = 3; % corresponds to the (symmetric) stable dist
initvec = [3 loc scale]; % [df loc scale]

ci_length_para = zeros([reps length(n_samp_vec)]);
coverage_para = zeros([reps length(n_samp_vec)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(['ex3: non-parametric bootstrap for tail_index = ', num2str(tail_index)]);
params = [tail_index, scale, loc];
[ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, dist, params, ES_stoy(n_tail_index), alpha, a);
disp(['Average nonparametric CI Length: ', num2str(average_length, ' % 7.4f')]);
disp(['Nonparametric Coverage Ratio:    ', num2str(mean_coverage_ratio, ' % 7.4f')]);

struct_nonpara_secondtailindex = struct('average_length', average_length, 'ci_length', ci_length, 'mean_coverage_ratio', mean_coverage_ratio, 'coverage_ratio', coverage_ratio);

%%%%%%%%%%%%%%%%%%%%%%%%
% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
ci_length = zeros([reps 1]);
coverage = zeros([reps 1]);

disp(delim); disp(['ex3: parametric bootstrap for tail_index = ', num2str(tail_index)]);
for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for i = 1:reps
        % generate random sample from a (symmetric) stable dist
        data = stabgen(n_samp_vec(k), tail_index, 0, scale ,loc); %don't set the seed, function takes random one
        initvec = [2 loc scale]; % [df loc scale]
        para_bs_hat = tlikmax(data, initvec);
        para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

        sim_para_mat = para_loc_hat + para_scale_hat * trnd(para_df_hat, [n_samp_vec(k) n_BS]);
        ES_vec = zeros([1 n_BS]);
        for j = 1:n_BS
            temp = sim_para_mat(:, j);
            params_temp = tlikmax(temp, initvec); %[df, loc, scale]

            % create bootstrap sample
            %ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
            %bs_samp = data(ind);

            % parametric
            %para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
            %para_loc_hat = para_bs_hat(1); para_scale_hat = para_bs_hat(2); para_df_hat = para_bs_hat(3);

            %para_bs_hat = tlikmax(bs_samp, initvec);
            %para_df_hat = para_bs_hat(1); para_loc_hat = para_bs_hat(2); para_scale_hat = para_bs_hat(3);

            % calculate theoretical ES based on the parameter estimates
            c01 = tinv(alpha, params_temp(1));
            ES_vec(j) = params_temp(2) + params_temp(3) * (-tpdf(c01, params_temp(1))/tcdf(c01, params_temp(1)) * (params_temp(1)+c01^2)/(params_temp(1)-1));
        end % j-loop

        % compute length of the CI and coverage
        ci_para = quantile(ES_vec, [a/2 1-a/2]);
        low_para = ci_para(1); high_para = ci_para(2);
        ci_length_para(i, k) = high_para - low_para;
        if ES_stoy(n_tail_index) >= low_para && ES_stoy(n_tail_index) <= high_para
            coverage_para(i, k) = 1;
        end

        if mod(i, 10) == 0
            disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
        end
    end % i-loop (reps)
end % k-loop (sample size)
disp(['Average parametric CI Length: ', num2str(mean(ci_length_para), ' % 7.4f')]);
disp(['parametric Coverage Ratio:    ', num2str(mean(coverage_para), ' % 7.4f')]);

%save 
struct_para_secondtailindex = struct('mean_ci_length_para', mean(ci_length_para), 'ci_length_para', ci_length_para, 'mean_coverage_ratio_para', mean(coverage_para), 'coverage_ratio_para', coverage_para);
struct_comb = struct('struct_nonpara_secondtailindex', struct_nonpara_secondtailindex, 'struct_para_secondtailindex', struct_para_secondtailindex);
save('ex3_secondtailindex_len+coverage.mat', 'struct_comb');

toc
disp(delim); disp(delim);  %took seconds
%% exercise 4 (calculate true ES of NCT (duplicate of ex2 part 2))
tic
% define parameters
delim = '************************************';
n_samp = 1e7; loc = 1; scale = 2;
%reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
df_vec = [3 6]; % degrees of freedom of the NCT
mu_vec = [-2 -1 0]; % (numerator) non-centrality parameter of the NCT
%theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
seed = rand*1000; alpha = .01; a = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% part 2 (calculate true ES of the NTC for the following parameters and via
%   (i) simulation and 
%   (ii) integral definition of the NTC for):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % two different df:                           |    3|    6|     |     |
% % four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|

ES_sim = zeros(numel(mu_vec), numel(df_vec));
ES_num = zeros(numel(mu_vec), numel(df_vec));

disp(delim); disp(['ES for different df and non-centrality', newline, 'parameters of the NCT (mu)']);
for df = 1:numel(df_vec)
        disp(delim); disp(['for df = ', num2str(df_vec(df))]);
    % (i) Simulation:
        for mu = 1:numel(mu_vec)
                ES_sim(mu, df) = loc + scale * Simulated_ES_NCT(n_samp, df_vec(df), mu_vec(mu), alpha, seed);
        end % end mu-loop
        disp(['via Simulation:          ', num2str(ES_sim(:, df)', '% 7.4f')]);

    % (ii) numeric integration:
        for mu = 1:numel(mu_vec)
                c01=nctinv(alpha , df_vec(df), mu_vec(mu));
                I01 = @(x) x.*nctpdf(x, df_vec(df), mu_vec(mu)); %note that the problem with nctpdf mentioned in footnote 11 on p.373 in the intermediate prob book has been solved in the standard matlab function, hence it is used here
                ES_num(mu, df) = loc + scale * integral(I01 , -Inf , c01) / alpha;
        end % end mu-loop
        disp(['via Numeric Integration: ', num2str(ES_num(:, df)', '% 7.4f')]);
end % df-loop (for df)

disp(delim);
struct_ES = struct('ES_sim', ES_sim, 'ES_num', ES_num);
save('ex2_trueES.mat', 'struct_ES');
% struct result:
%         | df = 3 | df = 6 |
% mu = -3 |        |        |
% mu = -2 |        |        |
% mu = -1 |        |        |
% mu = 0  |        |        |
toc
disp(delim); disp(delim);  %took 20.562306 seconds
%% exercise 4 (bootstrapping for first df)
% DGP = NCT and parametric bootstrap also NCT with own MLE funtion
% end outputs are the tables from Exercise 2
% two different df:                           |    3|    6|     |     |
% three different sample sizes (n_samp):      |  250|  500| 2000|     |
% four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|
tic
delim = '************************************';
loc = 1; scale = 2;
reps = 150; n_samp_vec = [250 2000]; n_BS = 1000; % note that n_samp = T
df = 3; % degrees of freedom of the NCT
n_df = 1;
mu_vec = [-2 -1 0]; % (numerator) non-centrality parameter of the NCT
%theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
%seed = rand*1000;
alpha = .01; a = 0.1; 

ci_length_nonpara = cat(3, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 3
coverage_ratio_nonpara = cat(3, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 

ci_length_para = cat(3, ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)]), ...
                        zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 3
coverage_para = cat(3, ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)]), ...
                       zeros([reps numel(n_samp_vec)])); % hardcode numel(mu_vec) = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(['for df = ', num2str(df)]);
for mu=1:numel(mu_vec)
    [ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 2, [scale loc df, mu_vec(mu)], ES_num(mu, n_df), alpha, a);
    disp(['Average nonparametric CI Length with mu = ', num2str(mu_vec(mu)),  ': ', num2str(average_length, '% 7.4f')]);
    disp(['Nonparametric Coverage Ratio: with mu =   ', num2str(mu_vec(mu)),  ': ', num2str(mean_coverage_ratio, '% 7.4f')]);
    ci_length_nonpara(:, :, mu) = ci_length;
    coverage_ratio_nonpara(:, :, mu) = coverage_ratio;
end % mu-loop

%%%%%%%%%%%%%%%%%%%%%%%%
%% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
ci_length = zeros([reps 1]);
coverage = zeros([reps 1]);

%!!check coverage rate (might be wrong still)!!
disp(delim); disp(['for df = ', num2str(df)]);
for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for mu = 1:numel(mu_vec)
    disp(['starting calculations for non-centrality param mu = ', num2str(mu_vec(mu))]);
        for i = 1:reps
            % generate random sample of a (regular) loc-scale t dist
            data = loc + scale * asymtrnd([n_samp_vec(k) 1], mu_vec(mu), df);
            initvec = [df mu_vec(mu) loc scale]; % [df loc scale]

            para_bs_hat = nctlikmax(data, initvec, 2); %the 2 is which pdf calculationmethod is taken
            para_df_hat = para_bs_hat(1); para_mu_hat = para_bs_hat(2); para_loc_hat = para_bs_hat(3); para_scale_hat = para_bs_hat(4);

            sim_para_mat = para_loc_hat + para_scale_hat * asymtrnd([n_samp_vec(k) n_BS], para_mu_hat, para_df_hat);
            ES_vec = zeros([1 n_BS]);
            for j = 1:n_BS
                %create bootstrap sample
                temp = sim_para_mat(:, j);
                %estimate parameters of bootstrap sample
                params_temp = nctlikmax(temp, initvec, 2); %[df mu loc scale]
                % calculate theoretical ES based on the parameter estimates
                c01 = nctinv(alpha , params_temp(1), params_temp(2));
                I01 = @(x) x.*nctpdf(x, params_temp(1), params_temp(2));
                ES_vec(j) = params_temp(3) + params_temp(4) * integral(I01 , -Inf , c01) / alpha;
            end % j-loop
            % compute length of the CI and coverage
            ci_para = quantile(ES_vec, [a/2 1-a/2]);
            low_para = ci_para(1); high_para = ci_para(2);
            ci_length_para(i, k, mu) = high_para - low_para;
            if ES_num(mu, n_df) >= low_para && ES_num(mu, n_df) <= high_para
                coverage_para(i, k, mu) = 1;
            end
            if mod(i, 10) == 0
                disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
            end
        end % i-loop (reps)
    end % mu-loop
    disp("mean length")
    disp(mean(ci_length_para));
    disp("mean coverage")
    disp(mean(coverage_para));
end % k-loop (samp size

% save
struct_nonpara_firstdf = struct('average_length', average_length, 'ci_length', ci_length, 'mean_coverage_ratio', mean_coverage_ratio, 'coverage_ratio', coverage_ratio);
struct_para_firstdf = struct('mean_ci_length_para', mean(ci_length_para), 'ci_length_para', ci_length_para, 'mean_coverage_ratio_para', mean(coverage_para), 'coverage_ratio_para', coverage_para);
struct_comb = struct('struct_nonpara_firstdf', struct_nonpara_firstdf, 'struct_para_firstdf', struct_para_firstdf);
save('results/ex4_firstdf_len+coverage.mat', 'struct_comb');

toc
disp(delim); disp(delim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exercise 4 (bootstrapping for second df = 6)
% DGP = NCT and parametric bootstrap also NCT with own MLE funtion
% end outputs are the tables from Exercise 2
% two different df:                           |    3|    6|     |     |
% three different sample sizes (n_samp):      |  250|  500| 2000|     |
% four different asymmetry parameters (mu):   |   -3|   -2|   -1|    0|
tic
delim = '************************************';
loc = 1; scale = 2;
reps = 150; n_samp_vec = [250 2000]; n_BS = 1000; % note that n_samp = T
df = 6; % degrees of freedom of the NCT
n_df = 2;
mu_vec = [-2 -1 0]; % (numerator) non-centrality parameter of the NCT
%theta = 0; % denominator non-centrality parameter of the NCT (for theta = 0 one gets the singly NCT)
%seed = rand*1000;
alpha = .01; a = 0.1; 

ci_length_nonpara = cat(3, ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]));
coverage_ratio_nonpara = cat(3, ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]));

ci_length_para = cat(3, ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]));
coverage_para = cat(3, ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]), ...
    zeros([reps numel(n_samp_vec)]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% non-parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(delim); disp(['for df = ', num2str(df)]);
for mu=1:numel(mu_vec)
    [ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 2, [scale loc df, mu_vec(mu)], ES_num(mu, n_df), alpha, a);
    disp(['Average nonparametric CI Length with mu = ', num2str(mu_vec(mu)),  ': ', num2str(average_length, '% 7.4f')]);
    disp(['Nonparametric Coverage Ratio: with mu =   ', num2str(mu_vec(mu)),  ': ', num2str(mean_coverage_ratio, '% 7.4f')]);
    ci_length_nonpara(:, :, mu) = ci_length;
    coverage_ratio_nonpara(:, :, mu) = coverage_ratio;
end % mu-loop

%%%%%%%%%%%%%%%%%%%%%%%%
%% parametric bootstrap %
%%%%%%%%%%%%%%%%%%%%%%%%
ci_length = zeros([reps 1]);
coverage = zeros([reps 1]);

%!!check coverage rate (might be wrong still)!!
disp(delim); disp(['for df = ', num2str(df)]);
for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    for mu = 1:numel(mu_vec)
    disp(['starting calculations for non-centrality param mu = ', num2str(mu_vec(mu))]);
        for i = 1:reps
            % generate random sample of a (regular) loc-scale t dist
            data = loc + scale * asymtrnd([n_samp_vec(k) 1], mu_vec(mu), df);
            initvec = [df mu_vec(mu) loc scale]; % [df loc scale]

            para_bs_hat = nctlikmax(data, initvec, 2); %the 2 is which pdf calculationmethod is taken
            para_df_hat = para_bs_hat(1); para_mu_hat = para_bs_hat(2); para_loc_hat = para_bs_hat(3); para_scale_hat = para_bs_hat(4);

            sim_para_mat = para_loc_hat + para_scale_hat * asymtrnd([n_samp_vec(k) n_BS], para_mu_hat, para_df_hat);
            ES_vec = zeros([1 n_BS]);
            for j = 1:n_BS
                %create bootstrap sample
                temp = sim_para_mat(:, j);
                %estimate parameters of bootstrap sample
                params_temp = nctlikmax(temp, initvec, 2); %[df mu loc scale]
                % calculate theoretical ES based on the parameter estimates
                c01 = nctinv(alpha , params_temp(1), params_temp(2));
                I01 = @(x) x.*nctpdf(x, params_temp(1), params_temp(2));
                ES_vec(j) = params_temp(3) + params_temp(4) * integral(I01 , -Inf , c01) / alpha;
            end % j-loop
            % compute length of the CI and coverage
            ci_para = quantile(ES_vec, [a/2 1-a/2]);
            low_para = ci_para(1); high_para = ci_para(2);
            ci_length_para(i, k, mu) = high_para - low_para;
            if ES_num(mu, n_df) >= low_para && ES_num(mu, n_df) <= high_para
                coverage_para(i, k, mu) = 1;
            end
            if mod(i, 10) == 0
                disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
            end
        end % i-loop (reps)
    end % mu-loop
    disp("mean length")
    disp(mean(ci_length_para));
    disp("mean coverage")
    disp(mean(coverage_para));
end % k-loop (samp size

% save
struct_nonpara_seconddf = struct('average_length', average_length, 'ci_length', ci_length, 'mean_coverage_ratio', mean_coverage_ratio, 'coverage_ratio', coverage_ratio);
struct_para_seconddf = struct('mean_ci_length_para', mean(ci_length_para), 'ci_length_para', ci_length_para, 'mean_coverage_ratio_para', mean(coverage_para), 'coverage_ratio_para', coverage_para);
struct_comb = struct('struct_nonpara_firstdf', struct_nonpara_seconddf, 'struct_para_firstdf', struct_para_seconddf);
save('results/ex4_firstdf_len+coverage.mat', 'struct_comb');

toc
disp(delim); disp(delim);