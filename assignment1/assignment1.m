%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Statistical Foundations for Finance - Homework 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question 1. compute pdf of the stable (first line: simulated density, second line: true density)
% set parameters of the pdf
a = 1.7; b = -0.4; c = 2; d = 0.3;
n = 40000; xvec = -20:.001:20;

%%% kernel density estimate %%%
% generate a random sample of size n from the S_{a,b}(c, d) distribution
% and plot the resulting density
randstab = stabgen(n, a, b, c, d, 2);
[f,x] = ksdensity(randstab, xvec);
figure, plot(x, f, 'g-', 'linewidth', 3)
xlim([-20 20])

%%% true density %%%
% calculate the actual theoretical values of a S_{a,b}(c, d) distribution
theostab = asymstabplus(xvec, a, b, c, d);
hold on, plot(xvec, theostab, 'b--', 'linewidth', 3), hold off

% prettyfy the plot by adding a legend, title, labels, etc
legend('Simulated PDF', ...
       'Theoretical PDF', ...
       'Location', 'southoutside', 'Orientation', 'vertical', 'Box', 'off');
title('PDFs For Stable Distribution')
xlabel("x"); ylabel("f(x)")
set(gca, 'fontsize', 10)

% save plot as .png file
saveas(gcf, 'assignment1_ex1.png')

%% Question 2. convolution of two independent stable random variables
% set parameters of the pdf
a = 1.7; b1 = -0.4; b2 = 1; c1 = 2; c2 = 1; d1 = -0.5; d2 = -0.3;
n = 40000; xvec = -20:.001:20;

% calculate the parameters of the stable r.v. that is given by the sum of two stable r.v.s
% (slide 535 in the lecture notes)
b_conv = (b1 * c1^a + b2 * c2^a)/(c1^a + c2^a); c_conv = (c1^a + c2^a)^(1/a); d_conv = d1 + d2;

%%% kernel density estimate %%%
% generate a random sample of size n from the S_{a, b}(d, c) distribution
% and plot the resulting density
randstab_conv = stabgen(n, a, b_conv, c_conv, d_conv, 2);
[f_conv,x_conv] = ksdensity(randstab_conv, xvec);
figure, plot(x_conv, f_conv, 'g-', 'linewidth', 3)
xlim([-20 20])

%%% true density %%%
% calculate the actual theoretical values of a S_{a, b}(c, d) distribution
theostab_conv = asymstabplus(xvec, a, b_conv, c_conv, d_conv);
hold on, plot(xvec, theostab_conv, 'b--', 'linewidth', 2), hold off

% prettyfy the plot by adding a legend, title, labels, etc
legend('Simulated PDF', ...
       'Theoretical PDF', 'Location', 'southoutside', 'Orientation', 'vertical', 'Box', 'off')
title('PDF For A Convolution of two Stable Distribution r.v.s')
xlabel("x"); ylabel("f(x)")
set(gca, 'fontsize', 10)

% save plot as .png file
saveas(gca, 'assignment1_ex2.png')

%% Question 3. convolution of two independent stable random variables with different tail index alpha
% set parameters of the pdf
% scale=c, location=d
%a1 = 1.6; a2 = 1.8;
a1 = 1.5; a2 = 1.9; % alphas for the second run
b = 0; c = 1; d = 0;
n = 4*1e2; xvec = -10:0.05:10; svec=-10:0.05:10;


%%% by simple integration formula: the convolution formula %%%
% slide "Asymmetric Stable: P.D.F. Calculation" (s. 553f)
f = convopdf(svec,a1,a2);
figure, plot(svec, f, 'g-', 'linewidth', 3)
xlim([-9.5 9.5])


%%% by the inversion formula applied to the characteristic function %%%
theostab_conv_ex3_2 = asymstabplus_ex3_2(xvec, a1, a2);
hold on, plot(xvec, theostab_conv_ex3_2, 'b--', 'linewidth', 3), hold off


%%% kernel density estimate %%%
n=10000;

% simulate a random sample of size n with different seeds to avoid
% correlation
randstab_conv_a1 = stabgen(n, a1, b, c, d, 5);
randstab_conv_a2 = stabgen(n, a2, b, c, d, 7585);
% add (convolute) the two random samples
randstab_conv_s = randstab_conv_a1 + randstab_conv_a2;
% apply the kernel density smoother in order to be able to plot the density
[f_conv_a, x_conv_a] = ksdensity(randstab_conv_s,svec);
hold on, plot(x_conv_a, f_conv_a, 'r:', 'linewidth', 3), hold off


%%% using the "conv" function in Matlab %%%
xvec = -10:.01:10;

% generate pdfs for the two alphas
s1 = asymstab_generalized(xvec, a1, b, c, d);
s2 = asymstab_generalized(xvec, a2, b, c, d);

% calculate the convolution
sconv = conv(s1, s2, 'same')/100;

% plot
hold on, plot(xvec', sconv, 'k:', 'LineWidth', 3), hold off

% prettyfy the plot by adding a legend, title, labels, etc
legend('PDF by convolution formula', ...
       'PDF by inversion formula', ...
       'PDF by simulation', ...
       'PDF by Matlab convolution formula', ...
       'Location', 'southoutside', 'Orientation', 'vertical', 'Box', 'off')
title('PDF For A Convolution of two Stable Distribution r.v.s with different \alpha')
xlabel("x"); ylabel("f(x)")
set(gca, 'fontsize', 10)

% save plot as .png file
if a1 == 1.6
    saveas(gca, 'assignment1_ex3.png')
else
    saveas(gca, 'assignment1_ex3_diffa.png')
end


%% Question 3(.2). convolution of THREE independent stable random variables with identical tail index alpha

a = 1.9; b1 = 0; b2 = -0.5; b3 = 0.4; c1 = 1; c2 = 1.5; c3 = 2; d1 = 0; d2 = -1; d3 = 1;
xvec = -16:.01:16; n=1e6; svec=-16:0.05:16;

% calculate the parameters of the stable r.v. that is given by the sum of two stable r.v.s
% (slide 535 in the lecture notes)
b_conv = (b1 * c1^a + b2 * c2^a + b3 * c3^a)/(c1^a + c2^a + c3^a);
c_conv = (c1^a + c2^a + c3^a)^(1/a);
d_conv = d1 + d2 + d3;


%%% kernel density estimate %%%
randstab_conv_s = stabgen(n, a, b_conv, c_conv, d_conv, 5);
% apply the kernel density smoother in order to be able to plot the density
[f_conv_a, x_conv_a] = ksdensity(randstab_conv_s, svec);
figure, plot(x_conv_a, f_conv_a, 'g-', 'linewidth', 3)
xlim([-15 15])

%%% using the "conv" function in Matlab %%%

% generate pdfs for the three alphas
s1 = asymstab_generalized(xvec, a, b1, c1, d1);
s2 = asymstab_generalized(xvec, a, b2, c2, d2);
s3 = asymstab_generalized(xvec, a, b3, c3, d3);

% calculate the convolution
sconv = conv(s3, conv(s1, s2, 'same')/100, 'same')/100;

% plot
hold on, plot(xvec', sconv, 'b--', 'LineWidth', 3), hold off

% prettyfy the plot by adding a legend, title, labels, etc
legend('PDF by simulation', ...
       'PDF by Matlab convolution formula', ...
       'Location', 'southoutside', 'Orientation', 'vertical', 'Box', 'off')
title('PDF For A Convolution of three Stable Distribution r.v.s with the same \alpha')
xlabel("x"); ylabel("f(x)")
set(gca, 'fontsize', 10)

% save plot as .png file
saveas(gca, 'assignment1_ex3_3rvs.png')

%% Question 4
a = 1.7; b = 0; c = 1; d = 0; xi = 0.01;

% theoretical ES using Stoyanov et al. (Book p. 490 - 492)
ES_stoy = asymstableES(xi, a, b, d, c, 1);
X = ['ES via Stoyanov et al: ', num2str(ES_stoy)]; 
disp(X);

% Simulation (Book p. 445)
nobs = 10^6;
ES_sim = Simulated_ES(nobs, a, b, c, d, xi, 0);
X = ['ES via simulation: ', num2str(ES_sim)]; 
disp(X);
% Result: -11.2002 with nobs = 10^6

% Evolution of ES 
% Loop over different observation sizes and compute the ES for each nobs
All_data = stabgen(10^9, a,b,c,d,0);
nobs = 10:1000:10^8;
ES_evolution = [];
for i=1:1000
    data = All_data(1:nobs(i));
    q = quantile(data, xi);
    Plo = data(data < q);
    ES_evolution(i) = mean(Plo);
end
nobs = 100000000:20000000:1000000000;
ES_evolution_2 = [];
for i=1:36
    data = All_data(1:nobs(i));
    q = quantile(data, xi);
    Plo = data(data < q);
    ES_evolution_2(i) = mean(Plo);
end

ES_total = [ES_evolution, ES_evolution_2];
nobs_total = [10:400000:10^8, 100000000:20000000:1000000000];
plot(nobs, ES_diff, 'linewidth', 2);
yline(ES_stoy, 'Linestyle', '--', 'Color', 'r', 'linewidth', 2);
xlabel('Number of Observations')
ylabel('Expected Shortfall')
ylim(-11.8, -11.1);
legend('ES via Simulation', 'location', 'NorthEast')

%% Question 5
a1 = 1.6; a2 = 1.8; %keep other parameters the same as before
xi = [0.01 0.025 0.05]; seed = 1; nobs = 1e6;

% Simulate the sum S = X1 + X2 and calculate ES for different values of xi

ES_sum_sim = Simulated_ES_sum(nobs, a1, a2, b, c, d, xi, seed);
X = ['ES via simulation for different levels: ', num2str(ES_sum_sim)]; 
disp(X);

% estimate parameters of the stable distribution of the sum. 
nobs = 1e6; X1 = stabgen(nobs, a1, b, c, d, 1); X2 = stabgen(nobs, a2, b, c, d, 2); S = X1 + X2;
[alpha,beta,sigma,mu] = stablereg(S);
X = ['Alpha: ', num2str(alpha), ' Beta: ', num2str(beta), ' Sigma: ', num2str(sigma), ' Mu: ', num2str(mu),]; 
disp(X);

% Now we can use the estimated parameters to calculate the ES via
% simulation
ES_stoy_sum = [];
for i = 1:3
    ES_stoy_sum(i) = asymstableES(xi(i), alpha, beta, mu, sigma ,1);
end
X = ['ES via Stoyanov et al: ', num2str(ES_stoy_sum)]; 
disp(X);

% Now we repeat it with new alphas
a1 = 1.5; a2 = 1.9; 
ES_sum_sim_2 = Simulated_ES_sum(nobs, a1, a2, b, c, d, xi, seed);
X = ['ES via simulation for different levels: ', num2str(ES_sum_sim_2)]; 
disp(X);

% Now using smaller set of simulated values (1e4) we estimate parameters of
% the stable distribution of the sum.
nobs = 1e6; X1 = stabgen(nobs, a1, b, c, d, 1); X2 = stabgen(nobs, a2, b, c, d, 2); S = X1 + X2;
[alpha,beta,sigma,mu] = stablereg(S);
X = ['Alpha: ', num2str(alpha), ' Beta: ', num2str(beta), ' Sigma: ', num2str(sigma), ' Mu: ', num2str(mu),]; 
disp(X);

% Now we can use the estimated parameters to calculate the 
ES_stoy_sum = [];
for i = 1:3
    ES_stoy_sum(i) = asymstableES(xi(i), alpha, beta, mu, sigma ,1);
end
X = ['ES via Stoyanov et al: ', num2str(ES_stoy_sum)]; 
disp(X);

