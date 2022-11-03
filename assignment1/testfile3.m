%testfile 3.1
a1=1.6;a2=1.8;
svec=-10:1:10;
f = convopdf(svec,a1,a2);
figure, plot(svec, f, 'linewidth', 3)

b=0;c=1;d=0;n=10000;
randstab_conv_a1 = stabgen(n, a1, b, c, d, 5);
randstab_conv_a2 = stabgen(n, a2, b, c, d, 7585);
randstab_conv_s = randstab_conv_a1 + randstab_conv_a2;
[f_conv_a, x_conv_a] = ksdensity(randstab_conv_s,svec);
hold on, plot(x_conv_a, f_conv_a, 'r--', 'linewidth', 3), hold off

legend('Theoretical PDF', ...
    'Simulated PDF', ...
    'Location', 'NorthWest');
title('PDFs For Sum of two Stable Distributions with different alphas')
%xlabel("x"); ylabel("S_{1.7, -0.4}(2, 0.3)(x)")
%set(gca, 'fontsize', 10)