%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just warm up by computing the ES of the location
%    scale Student t, and comparing analytic and
%    numeric integration.

% Set the tail prob alpha, and the Student t degrees
%    of freedom df, then calculate the left tail quantile
%    of the location 0, scale 1 distribution.
alpha=0.01; df=4; c01=tinv(alpha , df);

%% location zero, scale 1 ES for student t, calculating it
%   using first the analytic exact expression:
ES_01_analytic = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1) %#ok<*NOPTS>

% and now numeric integration:
I01 = @(x) x.*tpdf(x, df);
ES_01_numint = integral(I01 , -Inf , c01) / alpha

% they agree to about 14 digits!

%% now incorporate location and scale and check analytic vs numint
loc=1; scale=2; cLS = loc+scale*c01; % cLS is cutoff Location Scale
ES_wLS_analytic = loc+scale*ES_01_analytic
ILS = @(y) (y).*tpdf((y-loc)/scale, df)/scale;
ES_wLS_numint = integral(ILS , -Inf , cLS) / alpha
% again, perfect agreement.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inspect true ES versus simulation for a small sample size

% obtain the true ES, unknown in real life
loc=1; scale=2; df=4; alpha=0.01; % possibly new parameters
c01 = tinv(alpha , df); % left tail quantile, for loc-0 scale-1
truec = loc+scale*c01; % left tail quantile c
ES01 = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1);
trueES = loc+scale*ES01 % true theoretical ES

% simulate T-length vectors of IID location scale student t data
%   and record the empirical ES, and then compare to true
n_samp_vec=1e4; reps=200; ES_vec=zeros(reps,1);
for i=1:reps
  data=loc+scale*trnd(df,n_samp_vec,1); VaR=quantile (data, alpha);
  temp=data(data<=VaR); ES_vec(i)=mean(temp);
end

%% Now make a nice, visually appealing graphic:
figure
histogram(ES_vec), ax=axis;
set(gca,'fontsize', 8)
line ([ trueES trueES ] ,[0 ax(4)], 'color', 'g ', 'linewidth',3)
xlabel('ES value (simulation and true as vertical line)')
title(['Simulated Stud t Empirical ES, T=',int2str(n_samp_vec),' obs'])

%% exercise 1 (from first mail, 26.10.2022)
reps = 200; n_samp_vec = [250 500 2000]; n_BS = 1000; % note that n_samp = T
df = 2; loc = 1; scale = 2; alpha = .1;
initvec = [2 2 2]; % prior assumption for mle: [df, location, scale]

% set seed
rng(6, 'twister')

% initialize variables
ES_vec = zeros(n_BS, 1);
ci_length_para = zeros([reps length(n_samp_vec)]);
coverage_para = zeros([reps length(n_samp_vec)]);
% ci_length_nonpara = zeros([reps length(n_samp_vec)]);
% coverage_nonpara = zeros([reps length(n_samp_vec)]);

% calculate theoretical ES for the student t (analytically)
c01=tinv(alpha , df);
cLS = loc+scale*c01; % cLS is cutoff Location Scale
ES_01_analytic = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1);
ES_LS_analytic = loc+scale*ES_01_analytic;
trueES = ES_LS_analytic;

% nonparametric bootstrap
[average_length, coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 1, [scale, loc, df], trueES, alpha);
disp(['Average nonparametric CI Length: ', num2str(average_length)]);
disp(['Nonparametric Coverage Ratio: ', num2str(coverage_ratio)]);

% parametric bootstrap, choose which method for MLE should be used
mle_method = 1;
for k = 1:length(n_samp_vec)
    for i = 1:reps
        % generate "T" data points from the t-dist with "df" degrees of freedom
        % and location "loc" and scale "scale" (words in quotes refer to
        % variables)
        data = loc + scale * trnd(df, n_samp_vec(k), 1);
        
        % 2 methods for MLE (!! second methods gives wrong values!!)
        % i) matlab function mle
        if mle_method == 1
            para_bs_hat = mle(data, 'Distribution', 'tLocationScale'); % output: [loc scale df]
            para_loc_hat = para_bs_hat(1);
            para_scale_hat = para_bs_hat(2);
            para_df_hat = para_bs_hat(3);
        % ii) function from book 
        elseif mle_method == 2
            [param, stderr, iters, loglik, Varcov] = tlikmax(data, initvec);
            para_loc_hat = param(1);
            para_scale_hat = param(2);
            para_df_hat = param(3);
        end
        
        % generate parametric bootstrap sample with the estimated parameters
        
        % % generate parametric bootstrap sample with the estimated parameters
        %for j = 1:n_BS
        %   bs_samp = para_loc_hat + para_scale_hat * trnd(para_df_hat, n_samp_vec(k), 1);
        %   
        %   VaR = quantile(bs_samp, alpha/2);
        %   temp = bs_samp(bs_samp<=VaR);
        %   ES_vec(j) = mean(temp);
        %end
        %ci_para = quantile(ES_vec, [alpha/2 1-alpha/2]);
        %low_para = ci_para(1); high_para = ci_para(2);
        %ci_length_para(i, k) = high_para - low_para;
        %if ES_LS_analytic >= low_para && ES_LS_analytic <= high_para
        %    coverage_para(i,k) = 1;
        %end

        for j = 1:n_BS
            
            ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
            bs_samp = data(ind);
            
            % parametric
            para_bs_hat = mle(bs_samp, 'Distribution', 'tLocationScale'); %output: [loc scale df]
            para_bs_hat_MP = tlikmax(bs_samp, initvec) %MP: Marc Paolella
            %now calculate the theoretical ES based on the parameter
            %estimates
            c01=tinv(alpha , para_bs_hat_MP(1));
            ES_vec(j) = para_bs_hat_MP(2) + para_bs_hat_MP(3) * (-tpdf(c01,para_bs_hat_MP(1))/tcdf(c01,para_bs_hat_MP(1)) * (para_bs_hat_MP(1)+c01^2)/(para_bs_hat_MP(1)-1));
            disp(ES_vec(j))

            % non-parametric
            % directly compute "alpha"%-ES for each bootstrap sample
            % VaR = quantile(bs_samp, alpha/2);
            % temp=bs_samp(bs_samp<=VaR);
            % ES_vec(j)=mean(temp);
        end
        % compute length of the CI
        ci_para = quantile(ES_vec, [alpha/2 1-alpha/2]);
        low_para = ci_para(1); high_para = ci_para(2);
        ci_length_para(i, k) = high_para - low_para;
        if ES_LS_analytic >= low_para && ES_LS_analytic <= high_para
            coverage_para(i, k) = 1;
        end % j-loop (bootstrap loop)
    end %i-loop
end %k-loop

% compare the length of the CIs as a function of the sample size n_samp (= T)
%disp(['the average CI length with parametric bootstrap is:', num2str(mean(ci_length_para))])
%disp(['the average CI length with non-parametric bootstrap is:', num2str(mean(ci_length_nonpara))])
mean(ci_length_para)
% mean(ci_length_nonpara)

% compare coverage ratio of the theoretical ES as a function of the sample size n_samp (= T)
%disp(['the average coverage with parametric bootstrap is:', num2str(mean(coverage_para))])
%disp(['the average coverage with non-parametric bootstrap is:', num2str(mean(coverage_nonpara))])
mean(coverage_para)
% mean(coverage_nonpara)

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

% (i) non-parametric bootstrap
for df = 1:numel(df_vec)
    for mu=1:numel(mu_vec)
        [average_length, coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, 2, [df_vec(df), mu_vec(mu)], ES_num(mu), alpha);
        disp(['Average nonparametric CI Length with mu = ', num2str(mu_vec(mu)),  ': ', num2str(average_length)]);
        disp(['Nonparametric Coverage Ratio: with mu = ', num2str(mu_vec(mu)),  ': ' num2str(coverage_ratio)]);
    end % end mu
end

% (ii) parametric bootstrap based on student t distribution, !!Coverage rate is wrong!!
for df = 1:numel(df_vec)
    mle_method = 1;
    for k = 1:length(n_samp_vec)
        for mu = 1:numel(mu_vec)
            for i = 1:reps
                % generate "T" data points from the noncentral t-dist with "df" degrees of freedom
                % and location "loc" and scale "scale" and noncentrality parameter "mu" (words in quotes refer to
                % variables)
                data = loc + scale * asymtrnd(n_samp_vec(k), mu_vec(mu), df_vec(df), 1);
        
                % 2 methods for MLE (!! second methods gives wrong values!!)
                % We wrongly assume a t distribution and do MLE
                % (a) matlab function mle
                if mle_method == 1
                    para_bs_hat = mle(data, 'Distribution', 'tLocationScale'); % output: [loc scale df]
                    para_loc_hat = para_bs_hat(1);
                    para_scale_hat = para_bs_hat(2);
                    para_df_hat = para_bs_hat(3);
                % (b) function from book 
                elseif mle_method == 2
                    [param, stderr, iters, loglik, Varcov] = tlikmax(data, initvec);
                    para_loc_hat = param(1);
                    para_scale_hat = param(2);
                    para_df_hat = param(3);
                end
        
                % generate parametric bootstrap sample with the estimated parameters
        
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









