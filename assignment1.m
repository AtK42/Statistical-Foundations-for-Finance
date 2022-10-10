%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistical Foundations for Finance - Homework 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question 1. compute pdf of the stable (first line: simulated density, second line: true density)
a = 1.7; b = -0.4; c = 2; d = 0.3;
n = 40000; xvec = -20:.001:20;

%%% kernel density estimate %%%
% generate a random sample of size n from the S_{a,b}(d, c) distribution
% and plot the resulting density
randstab = stabgen(n, a, b, c, d, 2);
[f,x] = ksdensity(randstab, xvec);
figure, plot(x, f, 'r--', 'linewidth', 2)
xlim([-20 20])


%%% true density %%%
% calculate the actual theoretical values of a S_{a,b} distribution
theostab = asymstabplus(xvec, a, b, c, d);
hold on, plot(xvec, theostab, 'b-', 'linewidth', 2), hold off

pdfmatlab = makedist("Stable","alpha",a,"beta",b,"gam",c,"delta",d);
theostab1 = pdf(pdfmatlab, xvec);
hold on, plot(xvec, theostab1, 'c-', 'linewidth', 2), hold off

% prettyfy the plot
legend('Simulated PDF', ...
       'Theoretical PDF', ...
       'Theoretical PDF with Matlab', 'Location', 'NorthWest')
title('PDFs For Stable Distribution')
xlabel("x"); ylabel("S_{1.7, -0.4}(2, 0.3)(x)")
set(gca, 'fontsize', 10)

% 1. compute pdf of the stable (first line: kernel, second line: true density)
% define a sequence of x-values for which the pdf will be plotted and some
% different variables for the alpha parameter
%xvec = -10:.01:10;
%a1 = 0.8; a2 = 0.9; a3 = 1.1; a4 = 1.2; a5 = 1.9; b = 0.9;

% compute the pdf for different alpha parameters
%f1 = asymstab(xvec, a1, b); f2 = asymstab(xvec, a2, b); f3 = asymstab(xvec, a3, b); f4 = asymstab(xvec, a4, b); f5 = asymstab(xvec, a5, b);

% plot the pdfs
%legend_items = num2str([a1, a2, a3, a4, a5].', '\\alpha=%.1f');

%plot(xvec, f1, 'r-', xvec, f2, 'b-', xvec, f3, 'g-', xvec, f4 , 'r--', xvec, f5, 'b--', 'linewidth' , 2)
%legend ('\alpha = 0.8, \beta = 0.9', ...
%        '\alpha = 0.9, \beta = 0.9', ...
%        '\alpha = 1.1, \beta = 0.9', ...
%        '\alpha = 1.2, \beta = 0.9', ...
%        '\alpha = 1.9, \beta = 0.9', 'Location', 'NorthWest')

% simulate a stable random variable
% define variables
%nobs = 1e6;
%a = 1.2; b = 0.9; c = 0; d = 1;

% generate random sample of size nobs from the stable dist with the above
% defined variables
%randstab = stabgen(nobs, a, b, c, d, 42);

% plot the sample
%[f,x] = ksdensity(randstab);
%figure
%plot(x,f)

%% Question 2. convolution of two independent stable random variables
a = 1.7; b1 = -0.4; b2 = 1; c1 = 2; c2 = 2; d1 = -0.5; d2 = -0.3;
n = 1e5; xvec = -20:.001:20;

% by slide 535 in the lecture notes:
b_conv = (b1 * c1^a + b2 * c2^a)/(c1^a + c2^a); c_conv = (c1^a + c2^a)^(1/a); d_conv = d1 + d2;

%%% kernel density estimate %%%
% generate a random sample of size n from the S_{a,b}(d, c) distribution
% and plot the resulting density
randstab_conv = stabgen(n, a, b_conv, c_conv, d_conv, 2);
[f_conv,x_conv] = ksdensity(randstab_conv);
plot(x_conv, f_conv, 'r--', 'linewidth', 2)

%%% true density %%%
% calculate the actual theoretical values of a S_{a,b} distribution
theostab = asymstabplus(xvec, a, b_conv);
hold on, plot(xvec, theostab, 'b-', 'linewidth', 2), hold off

% prettyfy the plot
legend('Simulated PDF', ...
       'Theoretical PDF', 'Location', 'NorthWest')
title('PDFs For A Convolution of two Stable Distribution r.v.s')
xlim([-20 20])
xlabel("x"); ylabel("S(x)")
set(gca, 'fontsize', 10)

%% Question 4
a = 1.7; b = 0; c = 0; d = 1; xi = 0.01;

% theoretical ES using Stoyanov et al. (Book p. 490 - 492)
[ES_stoy, VaR] = asymstableES(xi, a, b, c, d,1);
%X = ['ES via Stoyanov et al: ', num2str(ES_stoy)]; 
%disp(X);


% Simulation (Book p. 445)
nobs = 10^6;
data = stabgen(nobs, a, b, d, c, 0);
q = quantile(data, xi);
Plo = data(data < q);
ES_sim = mean(Plo); 
X = ['ES via simulation: ', num2str(ES_sim)]; 
disp(X);
%Result: -10.981

%% Question 5

a1 = 1.6; a2 = 1.8; %keep other parameters the same as before
xi = [0.01 0.025 0.05]; seed = 0; nobs = 1e6;

% Simulate the sum S = X1 + X2 and calculate ES for different values of xi

ES_sum_sim = Simulated_ES(nobs, a1, a2, b, c, d, xi, seed);
X = ['ES via simulation for different levels: ', num2str(ES_sum_sim)]; 
disp(X);

% Now using smaller set of simulated values (1e4) we estimate parameters of
% the stable distribution of the sum.
nobs = 1e4; X1 = stabgen(nobs, a1, b, d, c, 0); X2 = stabgen(nobs, a2, b, d, c, 0); S = X1 + X2;
[alpha,beta,sigma,mu] = stablereg(S);
X = ['Alpha: ', num2str(alpha), ' Beta: ', num2str(beta), ' Sigma: ', num2str(sigma), ' Mu: ', num2str(mu),]; 
disp(X);

% Now we can use the estimated parameters to calculate the 
[ES_stoy, VaR] = asymstableES(xi, alpha, beta, mu, sigma ,1);
X = ['ES via Stoyanov et al: ', num2str(ES_stoy)]; 
disp(X);

% Now we repeat it with new alphas
a1 = 1.5; a2 = 1.9; 
ES_sum_sim = Simulated_ES(nobs, a1, a2, b, c, d, xi, seed);
X = ['ES via simulation for different levels: ', num2str(ES_sum_sim)]; 
disp(X);


