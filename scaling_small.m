%---------------------------------------------------------------------------------------------------
% For Paper
% "A Decomposition Approach to Multi-Agent Systems with Bernoulli Packet Loss"
% by Christian Hespe, Hamideh Saadabadi, Adwait Datar, Herbert Werner and Yang Tang
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

clear

addpath('graphs', 'performance')

%% Set up the problem
kappa = 0.1;
sysD  = ss(1, 1, 1, 0);
sysC  = ss(-kappa, 0, 0, 0);

p    = 0.5;
iter = 10;

% Iterate backwards to tackle the hard problems first. This improves the
% utilization of the CPU cores.
n = 20:-1:2;

%%
[N, ~] = meshgrid(n, 1:iter);

% Pre-allocate storage
norm_mean  = NaN(iter,length(n));
norm_enum  = NaN(iter,length(n));
norm_decom = NaN(iter,length(n));

stats_mean(iter,length(n))  = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);
stats_enum(iter,length(n))  = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);
stats_decom(iter,length(n)) = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);

tic
parfor i = 1:length(n)*iter
    L = full(laplace_matrix(circular_graph(N(i), 1, false)));
    
    [norm_mean(i),  ~, stats_mean(i) ] = performance_mean(sysD, sysC, L, p);
    [norm_enum(i),  ~, stats_enum(i) ] = performance_enumerated(sysD, sysC, L, p, true);
    [norm_decom(i), ~, stats_decom(i)] = performance_decomposed(sysD, sysC, L, p);
end
toc

% Calculate average timings
time_mean  = zeros(length(n),1);
time_enum  = zeros(length(n),1);
time_decom = zeros(length(n),1);
enum_prep  = zeros(length(n),1);
enum_solv  = zeros(length(n),1);
enum_yalm  = zeros(length(n),1);
for i = 1:length(n)
    time_mean(i)  = mean([stats_mean(:,i).total]);
    time_enum(i)  = mean([stats_enum(:,i).total]);
    time_decom(i) = mean([stats_decom(:,i).total]);
    enum_prep(i)  = mean([stats_enum(:,i).prep]);
    enum_solv(i)  = mean([stats_enum(:,i).solver]);
    enum_yalm(i)  = mean([stats_enum(:,i).yalmip]);
end

%% Plot data like in the paper
figure(1)
plot(n, norm_mean(1,:), 'k-.', n, norm_enum(1,:), 'k-', n, norm_decom(1,:), 'k--')
ylim padded
xlabel('Number of Agents')
ylabel('H_2-norm')
legend('Mean', 'Enumeration', 'Decomposed', 'Location', 'north west')

figure(2)
semilogy(n, time_mean, 'k-.', n, time_enum, 'k-', n, time_decom, 'k--')
ylim padded
xlabel('Number of Agents')
ylabel('Computation Time')
legend('Mean', 'Enumeration', 'Decomposed', 'Location', 'north west')

%% CSV Export
n = n';
norm_mean   = norm_mean(1,:)';
norm_enum   = norm_enum(1,:)';
norm_decom  = norm_decom(1,:)';
small_table = table(n, norm_mean, norm_enum, norm_decom, time_mean, time_enum, time_decom, enum_prep, enum_solv, enum_yalm);

name = sprintf('scaling_small_%d.csv', uint32(posixtime(datetime())));
writetable(small_table, name)
