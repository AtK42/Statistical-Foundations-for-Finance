function diff = stabcdfroot(x, xi, a, b)
    if exist('stableqkcdf.m','file')
        F = stableqkcdf(x, [a, b],1); % Nolan routine
    else
        [~, F] = asymstab1(x, a, b);
    end
    diff=F-xi;
end