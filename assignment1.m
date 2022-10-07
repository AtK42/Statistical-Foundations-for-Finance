%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistical Foundations for Finance - Homework 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exercise 1
%% compute pdf of the stable
% define a sequence of x-values for which the pdf will be plotted and some
% different variables for the alpha parameter
xvec = -10:.01:10;
a1 = 0.8; a2 = 0.9; a3 = 1.1; a4 = 1.2; a5 = 1.9; b = 0.9;

% compute the pdf for different alpha parameters
f1 = asymstab(xvec, a1, b); f2 = asymstab(xvec, a2, b); f3 = asymstab(xvec, a3, b); f4 = asymstab(xvec, a4, b); f5 = asymstab(xvec, a5, b);

% plot the pdfs
legend_items = num2str([a1, a2, a3, a4, a5].', '\\alpha=%.1f');
plot(xvec, f1, 'r-', xvec, f2, 'b-', xvec, f3, 'g-', xvec, f4 , 'y-', xvec, f5, 'b--', 'linewidth' , 1)
legend ('\alpha = 0.8, \beta = 0.9', ...
        '\alpha = 0.9, \beta = 0.9', ...
        '\alpha = 1.1, \beta = 0.9', ...
        '\alpha = 1.2, \beta = 0.9', ...
        '\alpha = 1.9, \beta = 0.9', 'Location', 'NorthWest')

%% simulate a stable random variable
% define variables
nobs = 1e6;
a = 1.2; b = 0.9; c = 0; d = 1;

% generate random sample of size nobs from the stable dist with the above
% defined variables
randstab = stabgen(nobs, a, b, c, d, 42);

% plot the sample
plot(1:1e6, randstab)
    
%% compute theoretial ES
a = 1.7; b = 0; c = 0; d = 1; xi = 0.01;


%% simulate ES
n = 10e6;
