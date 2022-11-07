function sample = asymtrnd(n_samp, loc, df)
% function to sample from the non-central t distribution
% returns an [n_samp 1] vector

% params
% % n_samp:         sample size
% % loc:            location parameter
% % df:             degrees of freedom

% Return NaN for elements corresponding to illegal parameter values.
% Non-integer degrees of freedom are allowed.
    df(df <= 0) = NaN;

    norm_rv = normrnd(loc, 1, [n_samp 1]);
    chisq_rv = chi2rnd(df, [n_samp 1]); % with 'normal' chisq dist
    noncent_chisq_rv = ncx2rnd(df, 0, [n_samp 1]); % with noncentral chisq dist
    
    sample = norm_rv / sqrt(chisq_rv/df) - norm_rv / sqrt(noncent_chisq_rv/df);

end % function