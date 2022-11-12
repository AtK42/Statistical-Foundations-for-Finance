function [ES, VaR] = asymstableES(xi , a, b, mu, scale , method)
    if nargin < 3, b=0; end, if nargin < 4, mu=0; end
    if nargin < 5, scale = 1; end, if nargin < 6, method = 1; end
    
    % Get q, the quantile from the S(0 ,1) distribution
    opt=optimset('Display', 'off', 'TolX', 1e-6);
    q=fzero(@stabcdfroot , -6, opt , xi , a, b);
    VaR=mu+scale*q;

    if(q==0)
        t0 = (1 / a) * atan(b * tan(pi * a/2));
        ES = ((2 * gamma((a-1) / a)) / (pi - 2*t0)) * (cos(t0) / cos(a*t0)^(1 / a));
        return ;
    end

    if(method==1)
        ES=(scale*Stoy(q,a,b)/xi)+mu;
    else
        ES=(scale*stabletailcomp (q,a,b) / xi )+ mu;
    end
end