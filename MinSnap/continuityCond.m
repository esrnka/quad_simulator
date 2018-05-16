function [Aeq_cont,beq_cont] = continuityCond(W)

% continuityCond : Returns Aeq, beq matrices that enforce continuity for a
%                  minimum snap problem
%
% INPUTS
%
% W ---------- Structure with the following elements:
%
%        tVecWp = Nx1 vector of time offsets from the initial time, in 
%                 seconds, with tVecWp(1) = 0. Each time offset represents
%                 the time to reach every waypoint
%
%             O = Polynomial order for each piecewise polynomial
%
% OUTPUTS
%
%   Aeq_cont = 4*(N-2)x(N-1)*(O+1) Aeq matrix that enforces continuity 
%              for a minimum snap problem
%
%   beq_cont = 4*(N-2)x1 beq matrix that enforces continuity for a
%              minimum snap problem
% 
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: 
%+==============================================================================+  

t_w = W.tVecWp;
n_w = numel(t_w);

%Observation: No boundary conditions on initial position and final pos
T0 = timeVectors(0,W.O);

%Allocate matrices
Aeq_cont = zeros(4*(n_w-2),(n_w-1)*(W.O+1));
beq_cont = zeros(4*(n_w-2),1);

for i = 2:n_w

    T = timeVectors(t_w(i) - t_w(i-1), W.O);
    
    %Enforce continuity
    if (i ~= n_w)
        range_row = 4*(i-2)+1:4*(i-1);
        range_col = (W.O + 1)*(i-2)+1:(W.O + 1)*i;
        Aeq_cont(range_row,range_col) = [T.Cp' -T0.Cp';
                                         T.Cv' -T0.Cv'; 
                                         T.Ca' -T0.Ca'; 
                                         T.Cj' -T0.Cj'];
        beq_cont(range_row,1) = zeros(4,1);
    end
end