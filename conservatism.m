clear

addpath('graphs', 'performance')

%% Set up the problem
kappa = 0.1;
sysD  = ss(1, 1, 1, 0);
sysC  = ss(-kappa, 0, 0, 0);

% Grid the probability axis. Higher density where larger changes are
% expected
p = [0, 1e-6, 1e-5, 1e-4, 1e-3:1e-3:9e-3, 0.01:0.01:0.29, 0.3:0.025:1];

% Calculate the norms for two examples, large and small networks
h = [4, 50];
N = h.*(h+1)/2;

%%
[H, P] = meshgrid(h, p);
sz = size(H);

% Pre-allocate storage
norm_mean  = NaN(length(p),length(h));
norm_enum  = NaN(length(p),length(h));
norm_decom = NaN(length(p),length(h));

tic
parfor i = 1:length(h)*length(p)
    L = full(laplace_matrix(triangle_graph(H(i))));
    
    norm_mean(i)  = performance_mean(sysD, sysC, L, P(i));
    norm_decom(i) = performance_decomposed(sysD, sysC, L, P(i));

    % The enumerated norm is only calculated for the small network, since
    % its scalling exponentially and thus infeasible to calculate for the
    % large network.
    [~, j] = ind2sub(sz, i);
    if j == 1
        norm_enum(i) = performance_enumerated(sysD, sysC, L, P(i), true);
    end
end
toc

%% Plot data like in the paper
figure(1)
plot(p, norm_mean(:,1), 'k-.', p, norm_enum(:,1), 'k-', p, norm_decom(:,1), 'k--')
ylim([0, 44])
title('Small Network Performance')
xlabel('Probability p')
ylabel('H_2-norm')
legend('Mean', 'Enumerated', 'Decomposed')

figure(2)
semilogy(p, norm_mean(:,2), 'k-.', p, norm_decom(:,2), 'k--')
ylim([3e1, 2e4])
title('Large Network Performance')
xlabel('Probability p')
ylabel('H_2-norm')
legend('Mean', 'Decomposed')

%% CSV Export
p = p';
small_mean  = norm_mean(:,1);
small_enum  = norm_enum(:,1);
small_decom = norm_decom(:,1);
small_table = table(p, small_mean, small_enum, small_decom);

large_mean  = norm_mean(:,2);
large_decom = norm_decom(:,2);
large_table = table(p, large_mean, large_decom);

name = sprintf('conservatism_small_%d.csv', uint32(posixtime(datetime())));
writetable(small_table, name)
name = sprintf('conservatism_large_%d.csv', uint32(posixtime(datetime())));
writetable(large_table, name)
