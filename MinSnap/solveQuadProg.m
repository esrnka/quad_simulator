function x = solveQuadProg(H, Aeq, beq)

% solveQuadProg : Solve a quadratic programming problem of the type
%                  min  x'Hx
%                  s.t. Aeq.x = beq
%                  
%                  Returns a vector x of length N_x
%
% INPUTS
%
%          H = N_x x N_x matrix of the minimization problem
%
%        Aeq = nr x N_x matrix for enforcing equality constraints
%
%        beq = nr x 1 vector for enforcing equality constraints
%
% OUTPUTS
%
%          x = N_x x 1 vector with optimal solution
% 
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: 
%+==============================================================================+
Nx = size(H,1);
Nr = size(beq,1);
%Form the M matrix
M = [H' Aeq'; Aeq, zeros(Nr,Nr)];
%Form the Y matrix
Y = [zeros(Nx,1); beq];
%Invert to solve
Z = M\Y;
%Extract the first Nx rows for the solution
x = Z(1:Nx,:);