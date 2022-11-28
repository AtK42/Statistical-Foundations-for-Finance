function pdfln = stdnctpdfln_j(x, nu, gam)
%pdf of NCT, computing the log density. x is the vector of points
%at which to evaluate the log density; nu and gam are the degrees of
%freedom and noncentrality parameters, respectively.
vn2 = (nu+1)/2; rho=x.^2;
pdfln = gammaln(vn2) - 1/2*log(pi*nu) - gammaln(nu/2) - vn2*log1p(rho/nu);
if(all(gam == 0)), return, end
idx = (pdfln >= -37); % -36.841 = log (1e-16)
if(any(idx))
    gcg = gam.^2; pdfln = pdfln - 0.5*gcg; xcg = x.*gam;
    term = 0.5*log(2) + log(xcg) - 0.5*log(max(realmin, nu+rho));
    term(term == -inf) = log(realmin); term(term == +inf) = log(realmax);
    maxiter = 1e4; k = 0;
    logterms = gammaln((nu+1+k)/2) - gammaln(k+1) - gammaln(vn2) + k*term;
    fractions = real(exp(logterms)); logsumk = log(fractions);
    while (k < maxiter)
        k = k + 1;
        logterms = gammaln((nu+1+k) / 2) - gammaln(k+1) - gammaln(vn2) + k*term(idx);
        fractions = real(exp(logterms-logsumk(idx)));
        logsumk(idx) = logsumk(idx) + log1p(fractions);
        idx(idx) = (abs(fractions) > 1e-4); if(all(idx == false)), break, end
    end
    pdfln = real(pdfln+logsumk);
end
%Program Listing 9.2: The direct density approximation to the NCT, computing the
%log density. x is the vector of points at which to evaluate the log density; nu and gam are
%the degrees of freedom and noncentrality parameters, respectively.