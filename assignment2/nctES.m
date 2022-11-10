function [ES, VaR] = nctES(xi, v, theta)
%code for quantile and ES given by the book (A.14, program listing A.2,
%page 467)
howfar = nctinv (1e-8, v, theta); % how far into the left tail to integrate
VaR = nctinv(xi, v, theta); % matlab routine for the quantile
I = quadl(@int, howfar, VaR,1e-6, [ ], v, theta); ES = I / xi;

function I = int(u, v, theta), pdf = nctpdf(u, v, theta); I = u.*pdf;