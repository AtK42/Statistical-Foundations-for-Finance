function [param,stderr,iters,loglik,Varcov] = weighted_MVTestimation3d(x, rho)
    % calculate weight vector
    T = length(x); tvec = (1:T); omega=(T-tvec + 1).^(rho-1); w = omega'/sum(omega);
    disp(w);
    % call MVTestimation function
    [param,stderr,iters,loglik,Varcov] = MVTestimation3d(x, w);
end