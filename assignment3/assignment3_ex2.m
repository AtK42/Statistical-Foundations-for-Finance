%% 2a
% only stuff to read

%% 2b
% implement function that calls MMF routine with hyperbolic weight decay
% "weighted_MMFAlgorithm"

%% 2c
% simulate a sequence of MVT r.v. 
T = 2000;
true_df = 4; % freely assumed
initial_df = .1; % freely assumed
dim = 3;
corr = corrmat(dim);
data = zeros([T,3]);
% make grid of T values that will be the df parameter
dfvec = linspace(6,3,T);

for t = 1:T
    % simulate 3d multivariate realizations
    data(t,:) = mvtrnd(corr, dfvec(t));
end

%% 2d
% create grid for rho values
reps = 1000
rhovec = linspace(0.2,1,9);
nu = zeros([1,length(rhovec)]);
mu = zeros([3,length(rhovec)]);

for i = 1:length(rhovec)
    [nu(:,i), nu_vec, mu(:,i), sigma] = weighted_MMFAlgorithm(rhovec(i), data, initial_df, reps);
end

% create plot !!TODO: Make plot nicer!!
plot(rhovec, nu)
title("Comparison rho versus estimated df")
xlabel("rho")
ylabel("estimated df")
save()

% as rho goes to 1 (equally weighted) the estimated df converges to 4.5
% which is the average df used.
% for rho = 0.2, estimated df is near 3 which makes sense because we give
% very high weight to the last few data points which have a df of 3 or a
% bit higher.

%%
% lower parameters than before
tic
reps = 200;
T = 200;
it = 200;
nu_mat = zeros([it, 9]);
for t = 1:it
    % generate new data
    corr = corrmat(dim);
    data = zeros([3,T]);
    for j = 1:T
        % simulate 3d multivariate realizations
        data(:,j) = mvtrnd(corr, dfvec(j));
    end
    for i = 1:length(rhovec)
        [nu_mat(t,i), nu_vec, mu, sigma] = weighted_MMFAlgorithm(rhovec(i), data, initial_df, reps);
    end
end
toc

%%
% create boxplots !!TODO: create nice graphic !!
C = mat2cell(nu_mat, [it], [1 1 1 1 1 1 1 1 1]);
figure('Position',[400 75 500 300])
boxplotGroup(C)
title('Boxplots for different values of rho', 'FontName','FixedWidth')

%% 2e
% repeat 2d but using brute-force MLE routine for the 3d MVT which needs to
% be extended to support the likelihood weight vector based on rho
tic
reps = 200;
T = 200;
it = 200;
nu_mat_bf = zeros([it, 9]);
for t = 1:it
    % generate new data
    corr = corrmat(dim);
    data = zeros([T,3]);
    for j = 1:T
        % simulate 3d multivariate realizations
        data(j,:) = mvtrnd(corr, dfvec(j));
    end
    for i = 1:length(rhovec)
        [param,stderr,iters,loglik,Varcov] = weighted_MVTestimation3d(data, rhovec(i));
        nu_mat_bf(t,i) = param(1);
    end
end
toc

%%
% create boxplots !!TODO: create nice graphic !!
C = mat2cell(nu_mat_bf, [20], [1 1 1 1 1 1 1 1 1]);
figure('Position',[400 75 500 300])
boxplotGroup(C)
title('Boxplots for different values of rho', 'FontName','FixedWidth')

%% 2f
% keep track of estimation time -> run everything in the end and note
% running times

