%---------------------------------------------------------------------------------------------------
% For Paper
% "A Decomposition Approach to Multi-Agent Systems with Bernoulli Packet Loss"
% by Christian Hespe, Hamideh Saadabadi, Adwait Datar, Herbert Werner and Yang Tang
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, Q, solver_stats] = performance_mean(sysD, sysC, L0, p)
%PERFORMANCE_MEAN Calculates a lower bound on the H2-norm of the given
%decomposable jump system with Bernoulli packet loss
%   This function uses Theorem 9 from the accompanying paper to calculate
%   a lower bound on the H2-norm of the given decomposable jump system.
%
%   Arguments:
%       sysD -> Decoupled part of the system
%       sysC -> Coupled part of the system
%       L0   -> Laplacian describing the nominal graph of the system
%       p    -> Probability of a successful package transmission
%   Returns:
%       H2           -> Lower bound on the H2-norm of the system
%       Q            -> Storage function matrix of the solution
%       solver_stats -> Timing information about the algorithm

% Setup the SDP solver
% The offset is used to convert strict LMI into nonstrict ones by pushing
% the solution away from the boundary.
offset = 1e-8;
opts   = sdpsettings('verbose', 0);

% Allocate storage for time measurements
solver_stats        = struct;
solver_stats.prep   = NaN;
solver_stats.solver = NaN;
solver_stats.total  = NaN;
solver_stats.yalmip = NaN;

% Initialize return values with reasonable defaults
H2    = inf;
timer = tic;

%% Prepare the system description
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ad, Bd, Cd, Dd] = ssdata(sysD);

nx = size(Ac,1);
nu = size(Bc,2);
N  = size(L0,1);

lambda = eig(L0);

%% Solve the H2-norm SDP as described in Theorem 9
Q    = sdpvar(nx, nx, N-1);
Z    = sdpvar(nu, nu, N-1);
cost = 0;
Constraints = [];

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = 2:N
    li = lambda(i);
    Qi = Q(:,:,i-1);
    Zi = Z(:,:,i-1);
    
    Acl = Ad + p*li*Ac;
    Bcl = Bd + p*li*Bc;
    Ccl = Cd + p*li*Cc;
    Dcl = Dd + p*li*Dc;
    
    LMI = Acl'*Qi*Acl + Ccl'*Ccl - Qi;
    TRC = Bcl'*Qi*Bcl + Dcl'*Dcl - Zi;
    cost = cost + trace(Zi);
    
    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        Constraints = [ Constraints                     ,...
                        Qi  >=  offset * eye(nx)        ,...
                        LMI <= -offset * eye(size(LMI)) ,...
                        TRC <= -offset * eye(size(TRC)) ];
    catch ME
        warning(ME.message)
        Q = [];
        return
    end
end
       
solver_stats.prep = toc(timer);
sol = optimize(Constraints, cost, opts);

if sol.problem ~= 0
    warning('YALMIP return an error: %s', sol.info)
    Q = [];
else
    % We calculate gamma^2 with the LMI constraints, so we need to take the
    % square root here.
    H2 = sqrt(value(cost));
    
    Q  = value(Q);
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
    solver_stats.total  = toc(timer);
end
end
