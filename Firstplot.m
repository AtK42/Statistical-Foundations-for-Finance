xvec = -6:0.02:6; b=-0.5;
a=1.8; f18=asymstab (xvec, a, b);% f28=asymstabpdf (xvec, a, b);
a=1.4; f14=asymstab (xvec, a, b);% f24=asymstabpdf (xvec, a, b);
figure, set(gca, 'fontsize',16)
plot (xvec, f14, '--r', xvec, f18, '--b', 'linewidth',1)
legend ( '\alpha = 1.4, \beta = -0.5', '\alpha = 1.8, \beta = -0.5', 'Location', 'NorthWest' )

%hold on, plot (xvec, f24 , ' g : ' , xvec , f28 , ' g : ' ) , hold of f , yl im ( [ 0 , 0 . 3 ] )
%figure , set (gca, 'fontsize',16)
%plot(xvec, f14-f24, 'k-'), title('PDF Differences for \alpha = 1.4')