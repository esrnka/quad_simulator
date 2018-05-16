function [Aeq_wp,beq_wp] = waypointCond(W)

% continuityCond : Returns Aeq, beq matrices that enforce waypoint
%                  nagivation in a minimum snap problem
%
% INPUTS
%
% W ---------- Structure with the following elements:
%
%        tVecWp = Nx1 vector of time offsets from the initial time, in 
%                 seconds, with tVecWp(1) = 0. Each time offset represents
%                 the time to reach every waypoint
%
%          rIwp = Nx3 matrix of waypoint CM positions in the I frame, in
%                 meters.  rIwp(k,:)' is the 3x1 waypoint position at time 
%                 tk = tVecWp(k).
%
%          vIwp = Nx3 matrix of waypoint CM velocities with respect to the
%                 I frame and expressed in the I frame, in meters/sec.
%                 vIwp(k,:)' is the 3x1 waypoint velocity at time
%                 tk = tVecWp(k). Set vIwp(k,:)' = NaN(3,1) if velocity is
%                 unconstrained through waypoint k
%
%          aIwp = Nx3 matrix of waypoint CM accelerations with respect to 
%                 the I frame and expressed in the I frame, in meters/sec^2.
%                 aIwp(k,:)' is the 3x1 waypoint acceleration at time tk =
%                 tVecWp(k). Set aIwp(k,:)' = NaN(3,1) if acceleration is
%                 unconstrained through waypoint k
%
%          jIwp = Nx3 matrix of waypoint CM jerk with respect to the I 
%                 frame and expressed in the I frame, in meters/sec^3.
%                 jIwp(k,:)' is the 3x1 waypoint jerk at time tk =
%                 tVecWp(k). Set jIwp(k,:)' = NaN(3,1) if jerk is
%                 unconstrained through waypoint k
%
%             O = Polynomial order for each piecewise polynomial
%
% OUTPUTS
%
% Aeq_wp ---------- Structure with the following elements:
%
%             x = 4*N x (N-1)*(O+1) Aeq matrix that enforces waypoint
%                 nagivation in the x direction in a minimum snap problem
%
%             y = 4*N x (N-1)*(O+1) Aeq matrix that enforces waypoint
%                 nagivation in the y direction in a minimum snap problem
%
%             z = 4*N x (N-1)*(O+1) Aeq matrix that enforces waypoint
%                 nagivation in the z direction in a minimum snap problem
%
% beq_wp ---------- Structure with the following elements:
%
%             x = 4*N x 1 beq matrix that enforces waypoint nagivation
%                 in the x direction in a minimum snap problem
%
%             y = 4*N x 1 beq matrix that enforces waypoint nagivation
%                 in the y direction in a minimum snap problem
%
%             z = 4*N x 1 beq matrix that enforces waypoint nagivation
%                 in the z direction in a minimum snap problem
% 
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: 
%+==============================================================================+ 

wp_x = [W.rIwp(:,1)'; W.vIwp(:,1)'; W.aIwp(:,1)'; W.jIwp(:,1)'];
wp_y = [W.rIwp(:,2)'; W.vIwp(:,2)'; W.aIwp(:,2)'; W.jIwp(:,2)'];
wp_z = [W.rIwp(:,3)'; W.vIwp(:,3)'; W.aIwp(:,3)'; W.jIwp(:,3)'];

t_w = W.tVecWp;
n_w = numel(t_w);
n = W.O;

%Observation: No boundary conditions on initial position and final pos
T0 = timeVectors(0,n);

%Allocate matrices
Aeq_wp.x = zeros(4*n_w,(n_w-1)*(n+1));
Aeq_wp.y = zeros(4*n_w,(n_w-1)*(n+1));
Aeq_wp.z = zeros(4*n_w,(n_w-1)*(n+1));

beq_wp.x = zeros(4*n_w,1);
beq_wp.y = zeros(4*n_w,1);
beq_wp.z = zeros(4*n_w,1);

%First waypoints
range_row = 1:4;
range_col = 1:n+1;
Aeq_wp.x(range_row,range_col) = [T0.Cp'; T0.Cv'; T0.Ca'; T0.Cj'];
Aeq_wp.y(range_row,range_col) = [T0.Cp'; T0.Cv'; T0.Ca'; T0.Cj'];
Aeq_wp.z(range_row,range_col) = [T0.Cp'; T0.Cv'; T0.Ca'; T0.Cj'];
beq_wp.x(range_row,1) = wp_x(:,1);
beq_wp.y(range_row,1) = wp_y(:,1);
beq_wp.z(range_row,1) = wp_z(:,1);

for i = 2:n_w

    T = timeVectors(t_w(i) - t_w(i-1),n);

    %Enforce waypoints
    range_row = 4*(i-1)+1:4*i;
    range_col = (n+1)*(i-2)+1:(n+1)*(i-1);
    Aeq_wp.x(range_row,range_col) = [T.Cp'; T.Cv'; T.Ca'; T.Cj'];
    Aeq_wp.y(range_row,range_col) = [T.Cp'; T.Cv'; T.Ca'; T.Cj'];
    Aeq_wp.z(range_row,range_col) = [T.Cp'; T.Cv'; T.Ca'; T.Cj'];
    beq_wp.x(range_row,1) = wp_x(:,i);
    beq_wp.y(range_row,1) = wp_y(:,i);
    beq_wp.z(range_row,1) = wp_z(:,i);
end

%Delete entries with NaN
nan_Index_x = find(isnan(beq_wp.x));
Aeq_wp.x(nan_Index_x,:) = [];
beq_wp.x(nan_Index_x) = [];

nan_Index_y = find(isnan(beq_wp.y));
Aeq_wp.y(nan_Index_y,:) = [];
beq_wp.y(nan_Index_y) = [];

nan_Index_z = find(isnan(beq_wp.z));
Aeq_wp.z(nan_Index_z,:) = [];
beq_wp.z(nan_Index_z) = [];
