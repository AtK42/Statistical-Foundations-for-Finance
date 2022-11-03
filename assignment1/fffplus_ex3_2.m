function I = fffplus_ex3_2(uvec, x, a1, a2, dopdf)
for ii=1:length(uvec)
    u=uvec(ii);
    t=(1-u)/u;
    cf=exp( - abs(t)^a1 - abs(t)^a2 );
    z=exp(-1i*t*x).*cf; if dopdf==1, g=real(z); else g=imag(z)./t;
    end
    I(ii)=g*u^(-2);
end