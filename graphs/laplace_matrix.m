%---------------------------------------------------------------------------------------------------
% For Paper
% "A Decomposition Approach to Multi-Agent Systems with Bernoulli Packet Loss"
% by Christian Hespe, Hamideh Saadabadi, Adwait Datar, Herbert Werner and Yang Tang
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function L = laplace_matrix(G)
%LAPLACE_MATRIX Calculate the Laplacian of the given graph
%   This function is a replacement for the Matlab built-in laplacian()
%   function that also works for digraphs, unlike the Matlab version.

A = adjacency(G);
D = diag(sum(A,2));
L = D - A;
end
