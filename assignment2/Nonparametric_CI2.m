function [ci_length, coverage_ratio, average_length, mean_coverage_ratio] = Nonparametric_CI2(reps, n_samp_vec, n_BS, dist, params, trueES, alpha, a)
% function to calculate the non-parametric CI

% input parameters:
% % reps:       number of repetitions
% % n_samp_vec: size of the random sample generated as a vector
% % n_BS:       number of bootstrap samples taken from each random sample
% % dist:       true distribution, can be either 
%               (i)     the symmetric student t (case 1)
%               (ii)    the asymmetric student t (case 2)
%               (iii)   the symmetric stable (case 3)
% % params      paremeter specifications for the true distribution and what the true ES would be
%               (i)     the symmetric student t:    param = [scale, location, df]
%               (ii)    the asymmetric student t:   param = [scale, location, df, mu]
%               (iii)   the symmetric stable:       param = [a, scale, location]

% output:       a vector with 
%               (i) lengths of the CIs (one for each repetition)
%               (ii) the coverage, i.e., percentage of the CIs that include the true ES

if dist == 1
    % check whether user input is valid
    if length(params) ~= 3 && isa(params, 'double') ~= 1
        error("'params' should be a double 1x3 vector where 'params(1)' is the scale, 'params(2)' is the location and 'params(3)' is the degrees of freedom")
    end
    % define each element of the params vector into a separate variable
    scale  = params(1); location = params(2); df = params(3);
elseif dist == 2
    % check whether user input is valid
    if length(params) ~= 2 && isa(params, 'double') ~= 1
        error("'params' should be a double 1x2 vector where 'params(1)' is the degrees of freedom and 'params(2)' is the (numerator) non-centrality parameter")
    end
    % define each element of the params vector into a separate variable
    scale = params(1); loc = params(2); df = params(3); mu = params(4);
elseif dist == 3
    % check whether user input is valid
    if length(params) ~= 3 && isa(params, 'double') ~= 1
        error("'params' should be a double 1x3 vector where 'params(1)' is the tail index alpha, 'params(2)' is the scale and 'params(3)' is the location")
    end
    % define each element of the params vecotr into a separate variable
    a_p = params(1); scale = params(2); location = params(3);
else
    error("Please specify a valid distribution. See function documentation for more information.");
end

% initialize variables
%average_length = zeros(length(n_samp),1);
%coverage_ratio = zeros(length(n_samp),1);
ES_vec = zeros(n_BS, 1);
ci_len = zeros([reps length(n_samp_vec)]);
coverage = zeros([reps length(n_samp_vec)]);

% set seed
rng(6, 'twister')

for k = 1:length(n_samp_vec)
    disp('***'); disp(['starting calculations for sample size = ', num2str(n_samp_vec(k))]);
    % initialize variables
    %len=zeros(reps,1);
    %coverage=zeros(reps,1);

    for i = 1:reps
    % first, generate the true dataset
        % (i) symmetric student t
        if dist == 1
            % reset seed (for reproducibility)
            % rng default;
            % generate the random sample
            data=location+scale*trnd(df,n_samp_vec(k),1); 
    
        % (ii) assymetric student t
        elseif dist == 2
            % generate the random sample
            data = loc + scale * asymtrnd(n_samp_vec(k), mu, df); % do not set a seed! Function creates random seeds

        % (iii) symmetric stable
        elseif dist == 3
            % generate the random sample
            data = stblrnd(a_p, 0, scale, location, n_samp_vec(k), 1);
        end

        %ESvec=zeros(n_BS,1);
        for j=1:n_BS
            % create bootstrap sample

            ind = unidrnd(n_samp_vec(k), [n_samp_vec(k) 1]);
            bs_samp = data(ind);

            % calculate ES
            VaR=quantile(bs_samp, alpha);
            temp=bs_samp(bs_samp<=VaR);
            ES_vec(j)=mean(temp);  
        end
        
        % calculate CI
        ci = quantile(ES_vec, [a/2 1-a/2]);
        lower_bound = ci(1); upper_bound = ci(2);
        ci_len(i, k) = upper_bound - lower_bound;

        % calculate Coverage
        if trueES >= lower_bound && trueES <= upper_bound
            coverage(i, k) = 1;
        end
        
        % info for user when running function
        if mod(i, 10) == 0
            disp(['finished rep ', num2str(i), ' out of ', num2str(reps), ' (' num2str(i/reps*100, '% 2.2f'), '% done)']);
        end
    end % i-loop (reps
end % end of k-loop

ci_length = ci_len;
coverage_ratio = coverage;
average_length = mean(ci_len);
mean_coverage_ratio = mean(coverage);

end % end of function


   