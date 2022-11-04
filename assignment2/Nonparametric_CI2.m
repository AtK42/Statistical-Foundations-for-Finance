function [length, coverage_ratio] = Nonparametric_CI2(reps, n_samp, n_BS, dist, params, trueES, alpha)
% function to calculate the non-parametric CI

% input parameters:
% % reps:       number of repetitions
% % n_samp:     size of the random sample generated
% % n_BS:       number of bootstrap samples taken from each random sample
% % dist:       true distribution, can be either 
%               (i)     the symmetric student t (case 1)
%               (ii)    the asymmetric student t (case 2)
%               (iii)   the symmetric stable (case 3)
% % params      paremeter specifications for the true distribution and what the true ES would be
%               (i)     the symmetric student t:    param = [scale, location, df]
%               (ii)    the asymmetric student t:   param = [df, mu]
%               (iii)   the symmetric stable:       param = [a, scale, location]

% output:       a vector with 
%               (i) lengths of the CIs (one for each repetition)
%               (ii) the coverage, i.e., percentage of the CIs that include the true ES


% initialize variables
length=zeros(reps,1);
coverage=zeros(reps,1);

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
    df = params(1); mu = params(2);
elseif dist == 3
    % check whether user input is valid
    if length(params) ~= 3 && isa(params, 'double') ~= 1
        error("'params' should be a double 1x3 vector where 'params(1)' is the tail index alpha, 'params(2)' is the scale and 'params(3)' is the location")
    end
    % define each element of the params vecotr into a separate variable
    a = params(1); scale = params(2); location = params(3);
else
    error("Please specify a valid distribution. See function documentation for more information.");
end

for i = 1:reps
    % first, generate the true dataset

    % (i) symmetric student t
    if dist == 1
        % reset seed (for reproducibility)
        rng default;
        % generate the random sample
        data=location+scale*trnd(df,n_samp,1); 
    
    % (ii) assymetric student t
    elseif dist == 2
        % generate the random sample
        data = asymt(df, mu, n_samp); % !!NEED TO CREATE THIS FUNCTION!!

    % (iii) symmetric stable
    elseif dist == 3
        % generate the random sample
        data = stblrnd(a, 0, scale, location, n_samp, 1);
    end

    ESvec=zeros(n_BS,1);
    for j=1:n_BS
      % create bootstrap sample
      bs_samp = datasample(data, n_samp);
      % calculate ES
      VaR=location+scale*quantile(bs_samp, alpha); temp=bs_samp(bs_samp<=VaR); ESvec(j)=mean(temp);  
    end

    % calculate CI
    lower_bound = quantile (ESvec, alpha); 
    upper_bound = quantile(ESvec, 1-alpha);
    length(i) = upper_bound - lower_bound;

    % calculate Coverage
    if trueES > lower_bound && trueES < upper_bound
        coverage(i) = 1;
    end

end % end of 'reps' loop

coverage_ratio = sum(coverage == 1)/reps;

end % end of function


   