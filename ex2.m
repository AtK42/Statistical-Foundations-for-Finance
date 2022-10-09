% exercise 2:
% formulas on page 282, but when alpha is equal, it is simple
% this is given
a = 1.7;
b1=-0.4;
c1=2;
d1=-0.5;

% we can choose b2,c2,d2
b2=-0.2;
c2=1.8;
d2=1;

% for the location it's just the sum
d=d1+d2;
d

% for the scale it's the sum of all c(i) to the power of alpha to the power
% of (1/alpha)
c=(c1^a + c2^a)^(1/a);
c

% for beta it's the sum of all betas times scales to the power of alpha over the
% the sum of scales to the power of alpha
b=(b1*c1^a + b2*c2^a)/(c1^a + c2^a);
b

f1= asymstabplus(-20:0.05:20,a,b1,c1,d1);
f2= asymstabplus(-20:0.05:20,a,b2,c2,d2);
f = asymstabplus(-20:0.05:20,a,b,c,d);

plot(-20:0.05:20, f,'DisplayName','f')
hold on
plot(-20:0.05:20, f1, 'DisplayName','f1')
plot(-20:0.05:20, f2,'DisplayName','f2')
hold off
legend


