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
figure, plot(x, f, 'r--', 'linewidth', 3)
xlim([-20 20])

%%% true density %%%
% calculate the actual theoretical values of a S_{a,b} distribution
theostab = asymstabplus(xvec, a, b, c, d);
hold on, plot(xvec, theostab, 'b-', 'linewidth', 3), hold off

%pdfmatlab = makedist("Stable","alpha",a,"beta",b,"gam",c,"delta",d);
%theostab1 = pdf(pdfmatlab, xvec);
%hold on, plot(xvec, theostab1, 'c-', 'linewidth', 2), hold off

% prettyfy the plot
legend('Simulated PDF', ...
       'Theoretical PDF', ...
       'Location', 'NorthWest');
title('PDFs For Stable Distribution')
xlabel("x"); ylabel("S_{1.7, -0.4}(2, 0.3)(x)")
set(gca, 'fontsize', 10)
saveas(gcf, 'assignment1_ex1.png')

%% Question 2. convolution of two independent stable random variables
a = 1.7; b1 = -0.4; b2 = 1; c1 = 2; c2 = 1; d1 = -0.5; d2 = -0.3;
n = 40000; xvec = -20:.001:20;

% by slide 535 in the lecture notes:
b_conv = (b1 * c1^a + b2 * c2^a)/(c1^a + c2^a); c_conv = (c1^a + c2^a)^(1/a); d_conv = d1 + d2;

%%% kernel density estimate %%%
% generate a random sample of size n from the S_{a,b}(d, c) distribution
% and plot the resulting density
randstab_conv = stabgen(n, a, b_conv, c_conv, d_conv, 2);
[f_conv,x_conv] = ksdensity(randstab_conv, xvec);
figure, plot(x_conv, f_conv, 'r--', 'linewidth', 3)
xlim([-20 20])
[f_conv,x_conv] = ksdensity(randstab_conv,xvec);
plot(x_conv, f_conv, 'r--', 'linewidth', 2)

%%% true density %%%
% calculate the actual theoretical values of a S_{a,b} distribution
theostab_conv = asymstabplus(xvec, a, b_conv, c_conv, d_conv);
hold on, plot(xvec, theostab_conv, 'b-', 'linewidth', 2), hold off

theostab = asymstabplus(xvec, a, b_conv, c_conv, d_conv);
hold on, plot(xvec, theostab, 'b-', 'linewidth', 2), hold off

% prettyfy the plot
legend('Simulated PDF', ...
       'Theoretical PDF', 'Location', 'NorthWest')
title('PDF For A Convolution of two Stable Distribution r.v.s')
xlabel("x"); ylabel("S(x)")
set(gca, 'fontsize', 10)
saveas(gca, 'assignment1_ex2.png')

%% Question 3. convolution of two independent stable random variables with different tail index alpha
% alpha is given as alpha1 = a1 = 1.6 and alpha2 = a2 = 1.8
% beta, scale and location are the same for both, i.e.
% beta=b=b1=b2=0 
% scale=c=c1=c2=1 
% location=d=d1=d2=0
a1 = 1.6; a2 = 1.8; b = 0; c = 1; d = 0;
n = 4*1e4; xvec = -10:.05:10;
svec=-10:0.05:10;


% now there are three different ways of computing the pdf for the
% convolution

% #1 Simple integration formula called the convolution formula, that I 
% showed you in class, along with formulae for difference, product, and
% quotient, remember? You program the integral convolution formula 
% (obviously, you need to be able to execute numeric integration), and 
% generate a plot of the resulting density (over, say, 400 points, or 
% however many you can do, i.e., maybe it takes too long, and you only do
% 100 points --- be smart, and understand what we are doing here.)
f = convopdf(svec,a1,a2);
figure, plot(svec, f, 'linewidth', 3)

%see slide "Asymmetric Stable: P.D.F. Calculation" (s. 553)
%figure, plot(-20:20, repelem(0, 41), 'g-', 'linewidth', 3) %place holder
%xlim([-20 20])

% #2 Next, you compute the pdf by using the inversion formula applied to 
% the characteristic function of the sum of X1 and X2, which is, remember,
% the product of the two characteristic functions. Do it, plot it. Overlay
% the two lines. They should be nearly identical.
theostab_conv_ex3_2 = asymstabplus_ex3_2(xvec, a1, a2);
hold on, plot(xvec, theostab_conv_ex3_2, 'b-', 'linewidth', 3), hold off

% #3 Add a third line to your graphic, based on simulation and kernel 
% density.
b=0;c=1;d=0;n=10000;
randstab_conv_a1 = stabgen(n, a1, b, c, d, 5);
randstab_conv_a2 = stabgen(n, a2, b, c, d, 7585);
randstab_conv_s = randstab_conv_a1 + randstab_conv_a2;
[f_conv_a, x_conv_a] = ksdensity(randstab_conv_s,svec);
hold on, plot(x_conv_a, f_conv_a, 'r--', 'linewidth', 3), hold off
 
% Obviously, lavishly annotate your graphic with titles, x and y labels,
% a legend, make the x and y numbers on the axis be big enough (use in
% Matlab set(gca,'fontsize',16)), etc. Use nice clear colors for each of
% the 3 lines of your graphic, say red, green, blue, or whatever you like
% (yellow is usually a bad choice), and also:

legend('PDF by convolution formula', ...
       'PDF by inversion formula', ...
       'PDF by simulation', ...
       'Location', 'NorthWest')
title('PDF For A Convolution of two Stable Distribution r.v.s with different \alpha')
xlabel("x"); ylabel("S(x)")
set(gca, 'fontsize', 10)
%saveas(gca, 'assignment1_ex3.png'

%% Question 4
a = 1.7; b = 0; c = 0; d = 1; xi = 0.01;

% theoretical ES using Stoyanov et al. (Book p. 490 - 492)
[ES_stoy, VaR] = asymstableES(xi, a, b, c, d,1);
X = ['ES via Stoyanov et al: ', num2str(ES_stoy)]; 
disp(X);


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


