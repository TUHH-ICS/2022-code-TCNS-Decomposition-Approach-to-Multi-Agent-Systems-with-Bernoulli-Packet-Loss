%---------------------------------------------------------------------------------------------------
% For Paper
% "A Decomposition Approach to Multi-Agent Systems with Bernoulli Packet Loss"
% by Christian Hespe, Hamideh Saadabadi, Adwait Datar, Herbert Werner and Yang Tang
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, Q, solver_stats] = performance_decomposed(sysD, sysC, L0, p)
%PERFORMANCE_DECOMPOSED Calculate an upper bound on the H2-norm of a
%decomposable jump system in a scalable manner
%   This function calculates an upper bound on the H2-norm of a
%   decomposable jump system in a scalable manner by applying a decoupled
%   analysis approach that is described in Theorem 7 of the accompanying
%   paper.
%
%   Arguments:
%       sysD -> Decoupled part of the system
%       sysC -> Coupled part of the system
%       L0   -> Laplacian describing the nominal graph of the system
%       p    -> Probability of a successful package transmission
%   Returns:
%       H2           -> Upper bound on the H2-norm of the system
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

%% Solve the SDP from Theorem 7
Q    = sdpvar(nx, nx);
Z    = sdpvar(nu, nu, N-1);
cost = 0;
Constraints = Q >= offset * eye(nx);

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = 2:N
    li  = lambda(i);
    lit = p*(1-p)*li;
    Zi  = Z(:,:,i-1);
    
    Acl = Ad + p*li*Ac;
    Bcl = Bd + p*li*Bc;
    Ccl = Cd + p*li*Cc;
    Dcl = Dd + p*li*Dc;
    
    LMI = Acl'*Q*Acl + Ccl'*Ccl - Q  + 2*lit * (Ac'*Q*Ac + Cc'*Cc);
    TRC = Bcl'*Q*Bcl + Dcl'*Dcl - Zi + 2*lit * (Bc'*Q*Bc + Dc'*Dc);
    cost = cost + trace(Zi);
    
    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        Constraints = [ Constraints                     ,...
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
