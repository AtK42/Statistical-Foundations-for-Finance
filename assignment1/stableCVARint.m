function [g] = stableCVARint(x, a, b)
if exist('stableqkpdf.m', 'file'), den = stableqkpdf(x,[a, b], 1);
else den = asymstab(x, a, b)';
end
g = x.*den;