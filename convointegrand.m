function f = convointegrand(xvec, s, a1, a2)
%function that computes the integrand used in the convolution formula
    f1 = asymstab(xvec, a1, 0);
    f2 = asymstab(s-xvec, a2, 0);
    f = (f1.*f2)';
end