%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistical Foundations for Finance - Homework 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. compute pdf of the stable (first line: simulated density, second line: true density)
a = 1.7; b = -0.4; c = 2; d = 0.3;
n = 1e4; xvec = -10:.001:10;

%%% kernel density estimate %%%
% generate a random sample of size n from the S_{a,b}(d, c) distribution
% and plot the resulting density
randstab = stabgen(n, a, b, c, d, 1);
[f,x] = ksdensity(randstab);
plot(x, f, 'r-', 'linewidth', 2)

%%% true density %%%
% calculate the actual theoretical values of a S_{a,b} distribution
theostab = asymstab(xvec, a, b);
hold on, plot(xvec, theostab, 'b-', 'linewidth', 2), hold off

%  prettyfy the plot
legend('Simulated PDF', ...
       'Theoretical PDF', 'Location', 'NorthWest')
title('PDFs For Stable Distribution')
xlim([-10 10])
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

%% 2.




%% Question 4
a = 1.7; b = 0; c = 0; d = 1; xi = 0.01;

% theoretical ES using Stoyanov et al. (Book p. 490 - 492)
%ES_stoy = asymstableES(xi, a, b, c, d,1);
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
%Result: -11.2002

%% Question 5
