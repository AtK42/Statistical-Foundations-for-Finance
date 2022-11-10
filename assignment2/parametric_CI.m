%alpha=0.1; df=4; loc=1; scale=2;
%T=250; reps=100; B=1000;
function [average_length, coverage_ratio] = parametric_CI(alpha, df, loc, scale, reps, T, B)
%True ES
trueES = ES_student_t(alpha,df,loc,scale);

%For parametric, this means you need to compute the MLE of the three 
% parameters corresponding to an IID student t DGP (thus: location, 
% scale, df) based on the bootstrap data set, and then compute the 
% theoretical ES based on those parameter values.

ci_len = zeros(reps,1);
coverage = zeros(reps, 1);

% simulate T-length vectors of IID location scale student t data
%   and record the empirical ES, and then compare to true
for i = 1:reps
    i
    data = loc + (scale * trnd(df,T,1));
    ind = unidrnd(T, [T, B]);
    ES_vec = zeros(B,1);
    for k = 1:B
        bs_samp = data(ind(:,k));
        mle_t = tlikmax0(bs_samp,[df loc scale]);
        ES_vec(k) = ES_student_t(alpha, mle_t(1), mle_t(2), mle_t(3));
    end
    % calculate CI
    ci = quantile(ES_vec, [alpha/2 1-alpha/2]);
    lower_bound = ci(1); upper_bound = ci(2);
    ci_len(i) = upper_bound - lower_bound;
    % calculate Coverage
    if trueES >= lower_bound && trueES <= upper_bound
        coverage(i) = 1;
    end
end
format long
average_length = mean(ci_len);
coverage_ratio = mean(coverage);

end