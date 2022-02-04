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
h = 142:-2:2;
N = h.*(h+1)/2;

%%
[H, ~] = meshgrid(h, 1:iter);

% Pre-allocate storage
norm_mean  = NaN(iter,length(h));
norm_decom = NaN(iter,length(h));

stats_mean(iter,length(h))  = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);
stats_decom(iter,length(h)) = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);

tic
parfor i = 1:length(h)*iter
    L = full(laplace_matrix(triangle_graph(H(i))));
    
    [norm_mean(i),  ~, stats_mean(i) ] = performance_mean(sysD, sysC, L, p);
    [norm_decom(i), ~, stats_decom(i)] = performance_decomposed(sysD, sysC, L, p);
end
toc

% Calculate average timings
time_mean  = zeros(length(h),1);
time_decom = zeros(length(h),1);
for i = 1:length(h)
    time_mean(i)  = mean([stats_mean(:,i).total]);
    time_decom(i) = mean([stats_decom(:,i).total]);
end

%% Plot data like in the paper
figure(1)
loglog(N, norm_mean(1,:), 'k-.', N, norm_decom(1,:), 'k--')
ylim padded
xlabel('Number of Agents')
ylabel('H_2-norm')
legend('Mean', 'Decomposed', 'Location', 'north west')

figure(2)
loglog(N, time_mean, 'k-.', N, time_decom, 'k--')
ylim padded
xlabel('Number of Agents')
ylabel('Computation Time')
legend('Mean', 'Decomposed', 'Location', 'north west')

%% CSV Export
h = h';
N = N';
norm_mean  = norm_mean(1,:)';
norm_decom = norm_decom(1,:)';
large_table = table(h, N, norm_mean, norm_decom, time_mean, time_decom);

name = sprintf('scaling_large_%d.csv', uint32(posixtime(datetime())));
writetable(large_table, name)
