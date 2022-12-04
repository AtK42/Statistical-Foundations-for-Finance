function [v,mu, sigma] = Algo(X,epsilon, w)
%Modified EM algorithm thing
 %T=10000; C = randcorr(3); df = 4;
 %x = mvtrnd(C, df, T)';
 %w = ones(length(x) , 1)*(1/length(x));
 %epsilon = 0.05;

d = size(X, 1);
v = epsilon;
mu = 1/length(X)*sum(X, 2);
%initialize sigma matrix
sigma0 = 0;
for i = 1:length(X)
    sigma0 = sigma0 + (X(:, i) - mu)*(X(:, i) - mu)';
end


sigma = sigma0/length(X);
%gammaold = 100;
%sigmaold = 100;
%muold = 100;
%gamma = 200;
r = 1;
while (r < 500)
    %E-step, compute weights
    %gammaold = gamma;
    %muold = mu;
    %sigmaold = sigma;
    [delta gamma] = calculateGamma(X, v, d, sigma, mu);
    
    mu = updateMu(w, gamma, X);
    sigma = updateSigma(w, gamma, X, mu);
    v = updateV(w, X, mu, v, d, sigma);
    r = r + 1
end

end
function [mu] = updateMu(w, gamma, X)
numerator = 0;
denominator = 0;
for i = 1:length(X)
    numerator = numerator + w(i)*gamma(i)*X(:, i);
    denominator = denominator + w(i)*gamma(i);
end
mu = numerator/denominator;
end

function [sigma] = updateSigma(w, gamma, X, mu)
sigma = 0;
denominatorSum = 0;
    for j = 1:length(X)
        denominatorSum = denominatorSum + w(j)*gamma(j);
    end
    for i = 1:length(X)
        sigma = sigma + (w(i)*gamma(i)*(X(:,i)-mu)*(X(:,i)-mu)')/denominatorSum;
    end
end

function [v] = updateV(w, X, mu, oldV, d, sigma)
 vSum = 0;
 n = length(X);
 [delta , gamma] = calculateGamma(X, oldV, d, sigma, mu);
    for i = 1:length(X)
    vSum = vSum + w(i) * ((oldV + d)/(oldV + delta(i)) - log((oldV + d)/(oldV + delta(i)))- 1);
    end
 v = fzero(@(x)  calculatePhi(x/2) - calculatePhi((x+d)/2) + vSum , [10^-6, 10^15]);
end

function [phi] = calculatePhi(x)
phi = psi(x) - log(x);
end

function [delta gamma] = calculateGamma(X, v, d, sigma, mu)
delta = zeros(length(X), 1);
gamma = zeros(length(X), 1);
    for i = 1:length(X)
       delta(i, 1) = (X(:, i) - mu)' * inv(sigma) * (X(:, i) - mu);
       gamma(i, 1) = (v + d)/(v + delta(i,1));
    end
end