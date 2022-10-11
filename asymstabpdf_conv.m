function f = asymstabpdf_conv (xvec, yvec, a1, a2, b)
% attempt to get a numerical solution to the sum of two stable r.v.s

% for a1 or a2 == 1
if a1==1, error('not ready yet'), end
if a2==1, error('not ready yet'), end

% define theta0
t0 = (1/a)*atan(b*tan((pi*a)/2));
if x == 0; gamma(1+1/a) * cos(t0) / ( pi*(1+ (b*tan(pi*a/2))^2)^(1/(2*a)) ); end
