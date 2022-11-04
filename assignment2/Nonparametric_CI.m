
function [Length, Coverage_ratio] = Nonparametric_CI(T, rep, B, dist, param, trueES, alpha)
% This function takes parameters for the number of simulations, true distribution (see below), the 
% parameters for the true distribution and what the true ES would be. 

% The output is a vector with all lengths of the CI's and the coverage, i.e., how many percent of the
% CI's include the true ES

% parameter specifications:
% student-t (case 1): param = [scale, location, df]
% singly noncentral student-t (case 2): param = [k, mu]
% symmetric stable (case 3): param = [a, scale, location]


Length=zeros(rep,1);
Coverage=zeros(rep,1);

for i = 1:rep

    % generate true dataset
    if dist == 1
        % symmetric student t
        scale  = param(1); location = param(2); df = param(3);
        data=location+scale*trnd(df,T,1); 
    
    elseif dist == 2
        % assymetric student t
        k = param(1); mu = param(2);
        data = asymt(k, mu, T); % !!NEED TO CREATE THIS FUNCTION!!

    elseif dist == 3
        % symmetric stable paretian
        a = param(1); scale = param(2); location = param(3);
        data = stblrnd(a, 0, scale, location, T, 1);
    else
        error("Need to specify a distribution!");
    end

    ESvec=zeros(B,1);
    for j=1:B
      % create bootstrap sample
      bootstrap_sample = datasample(data, T);
      % calculate ES
      VaR=location+scale*quantile(bootstrap_sample, alpha); temp=bootstrap_sample(bootstrap_sample<=VaR); ESvec(j)=mean(temp);  
    end

    % calculate CI
    lower_bound = quantile (ESvec, alpha); 
    upper_bound = quantile(ESvec, 1-alpha);
    Length(i) = upper_bound - lower_bound;

    % calculate Coverage
    if trueES > lower_bound && trueES < upper_bound
        Coverage(i) = 1;
    end
end
Coverage_ratio = sum(Coverage == 1)/rep;
end



   