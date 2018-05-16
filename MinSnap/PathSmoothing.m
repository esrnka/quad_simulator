function [R] = PathSmoothing(W)
% PathSmoothing_temp : Generates a smooth trajectory from a set of 
%                      waypoints.
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
%            Ts = Desired sampling frequency for the output trajectory
%
% OUTPUTS
%
% R ---------- Structure with the following elements:
%
%          tVec = Nx1 vector of uniformly-sampled time offsets from the
%                 initial time, in seconds, with tVec(1) = 0.
%
%        rIstar = Nx3 matrix of desired CM positions in the I frame, in
%                 meters.  rIstar(k,:)' is the 3x1 position at time tk =
%                 tVec(k).
%
%        vIstar = Nx3 matrix of desired CM velocities with respect to the I
%                 frame and expressed in the I frame, in meters/sec.
%                 vIstar(k,:)' is the 3x1 velocity at time tk = tVec(k).
%
%        aIstar = Nx3 matrix of desired CM accelerations with respect to the I
%                 frame and expressed in the I frame, in meters/sec^2.
%                 aIstar(k,:)' is the 3x1 acceleration at time tk =
%                 tVec(k).
%
%        jIstar = Nx3 matrix of desired CM jerk with respect to the I
%                 frame and expressed in the I frame, in meters/sec^3.
%                 jIstar(k,:)' is the 3x1 jerk at time tk = tVec(k).
%
%       optCost = 1x1 constant evaluating the optimal minimum snap cost
%
%+------------------------------------------------------------------------------+
% References: http://ieeexplore.ieee.org/document/5980409/
%
%
% Author: Marcelino Mendes de Almeida
%+==============================================================================+  

%Quadratic matrix
H = quadraticMatrix(W);

%Continuity enforcement
[Aeq_cont,beq_cont] = continuityCond(W);

%Enforce waypoints
[Aeq_wp,beq_wp] = waypointCond(W);

%% Build optimization problem
Aeq_x = [Aeq_cont; Aeq_wp.x];
Aeq_y = [Aeq_cont; Aeq_wp.y];
Aeq_z = [Aeq_cont; Aeq_wp.z];
beq_x = [beq_cont; beq_wp.x];
beq_y = [beq_cont; beq_wp.y];
beq_z = [beq_cont; beq_wp.z];

%Solve for x
sol_ax = solveQuadProg(H, Aeq_x, beq_x);

%Solve for y
sol_ay = solveQuadProg(H, Aeq_y, beq_y);

%Solve for z
sol_az = solveQuadProg(H, Aeq_z, beq_z);

%Get the optimal cost
W.optCost = sol_ax'*H*sol_ax + sol_ay'*H*sol_ay + sol_az'*H*sol_az;

%% Get final values of pvaj
n_w = numel(W.tVecWp);  % Number of waypoints
sol_ax = reshape(sol_ax,[W.O+1 (n_w-1)]);
sol_ay = reshape(sol_ay,[W.O+1 (n_w-1)]);
sol_az = reshape(sol_az,[W.O+1 (n_w-1)]);
dt = W.Ts;

R.tVec = W.tVecWp(1):dt:W.tVecWp(end);
k = 1;
n_pts = numel(R.tVec);
R.rIstar = zeros(n_pts, 3);
R.vIstar = zeros(n_pts, 3);
R.aIstar = zeros(n_pts, 3);
R.jIstar = zeros(n_pts, 3);
for i = 1:n_pts
    if (R.tVec(i) > W.tVecWp(k+1))
        k = k+1;
    end

    T = timeVectors(R.tVec(i) - W.tVecWp(k), W.O);
    R.rIstar(i,:) = T.Cp'*[sol_ax(:,k) sol_ay(:,k) sol_az(:,k)];
    R.vIstar(i,:) = T.Cv'*[sol_ax(:,k) sol_ay(:,k) sol_az(:,k)];
    R.aIstar(i,:) = T.Ca'*[sol_ax(:,k) sol_ay(:,k) sol_az(:,k)];
    R.jIstar(i,:) = T.Cj'*[sol_ax(:,k) sol_ay(:,k) sol_az(:,k)];
end