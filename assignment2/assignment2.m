%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Just warm up by computing the ES of the location
%    scale Student t, and comparing analytic and
%    numeric integration.

% Set the tail prob alpha, and the Student t degrees
%    of freedom df, then calculate the left tail quantile
%    of the location 0, scale 1 distribution.
alpha=0.01; df=4; c01=tinv(alpha , df);

%% location zero, scale 1 ES for student t, calculating it
%   using first the analytic exact expression:
ES_01_analytic = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1) %#ok<*NOPTS>
% and now numeric integration:
I01 = @(x) x.*tpdf(x, df);
ES_01_numint = integral(I01 , -Inf , c01) / alpha
% they agree to about 14 digits!

%% now incorporate location and scale and check analytic vs numint
loc=1; scale=2; cLS = loc+scale*c01; % cLS is cutoff Location Scale
ES_wLS_analytic = loc+scale*ES_01_analytic
ILS = @(y) (y).*tpdf((y-loc)/scale, df)/scale;
ES_wLS_numint = integral(ILS , -Inf , cLS) / alpha
% again, perfect agreement.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inspect true ES versus simulation for a small sample size

% obtain the true ES, unknown in real life
loc=1; scale=2; df=4; alpha=0.01; % possibly new parameters
c01 = tinv(alpha , df); % left tail quantile, for loc-0 scale-1
truec = loc+scale*c01; % left tail quantile c
ES01 = -tpdf(c01,df)/tcdf(c01,df) * (df+c01^2)/(df-1);
trueES = loc+scale*ES01 % true theoretical ES

% simulate T-length vectors of IID location scale student t data
%   and record the empirical ES, and then compare to true
T=1e4; rep=200; ESvec=zeros(rep,1);
for i=1:rep
  data=loc+scale*trnd(df,T,1); VaR=quantile (data, alpha);
  temp=data(data<=VaR); ESvec(i)=mean(temp);
end

%% Now make a nice, visually appealing graphic:
figure
histogram(ESvec), ax=axis;
set(gca,'fontsize', 8)
line ([ trueES trueES ] ,[0 ax(4)], 'color', 'g ', 'linewidth',3)
xlabel('ES value (simulation and true as vertical line)')
title(['Simulated Stud t Empirical ES, T=',int2str(T),' obs'])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% have a look at makedist
%For the data generating process (DGP) of IID location-scale
%  Student t, you fix the 3 true parameters, like we did above,
%  and then in a FOR loop, you simulate "rep" repetitions of
%  an IID T-length sequence of Student t, say rep=1000.
%  (Naturally, you try your codes with rep=2, to ensure they
%  work, and then choose rep based on your computing resources, e.g.,
%  perhaps 1000 is too high for your little laptop...)

%  Make your program general in that T and rep are either defined
%  up front, or passed as an argument.

%Now, for each of those "rep" repetitions, you calculate a bootstrap
%  confidence interval, 90%, based on B bootstrap replications, where
%  B should be at least 1000, though you can try smaller values,
%  and assess if, say, B=250 is in fact enough.

%For each of the "rep" CIs, you then compute two things:
%    a) the length of the interval,
%    b) whether or not the interval contains the TRUE ES, which you
%        of course have, because we know the true DGP. The above
%        codes have the computation of the theoretical Student t
%       ES in there, as you saw.  Notice the vector produced is a
%       sequence of zeros and ones (or "true" and "false", for you
%       computer science people...).

rep = 10; T = 1e4; B = 250;
df = 2; loc = 1; scale = 2; alpha = .1;

range default
for i = 1:rep
    data = loc + scale * trnd(df, T, 1);
    ind = unidrnd(T, [T B]); bs_samp = data(ind);
    ci = quantile(bs_samp, [alpha/2 1-alpha/2]); low = ci(1); high = ci(2);
    ci_length = high - low; 
end