function [H] = quadraticMatrix(W)

% quadraticMatrix : Returns the H matrix of the cost function for a minimum
%                   snap problem
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
%             H = (N-1)*(O+1) x (N-1)*(O+1) matrix for solving the
%                 quadratic programming problem
%
%+------------------------------------------------------------------------------+
% Author: Evan Srnka
%+==============================================================================+
%Determine Hi matrix size
N = length(W.tVecWp);
Hi = zeros(W.O+1,W.O+1);

%Determine delta t vector
delt = [];
for h = 1:N-1
    delt = [delt, W.tVecWp(h+1)-W.tVecWp(h)];
end

%Populate the Hi matrix
for i = 1:N-1
    Hi(:,:,i) = zeros(W.O+1,W.O+1);
    for j = 0:(W.O)
        for k = 0:(W.O)
            if j <= 3 || k <= 3
                Hi(j+1,k+1,i) = 0;
            else
                alpha = j+k-7;
                B = j*(j-1)*(j-2)*(j-3)*k*(k-1)*(k-2)*(k-3)*(1/alpha);
                Hi(j+1,k+1,i) = B*delt(i)^(alpha);
            end
        end
    end
end

%Build the big H matrix
H = Hi(:,:,1);
for l = 2:size(Hi,3)
    H = blkdiag(H,Hi(:,:,l));
end

                             
end