function [R] = pathplan(delt, Tsim, plotFlag)
% stateEstimatorUKF : Unscented-Kalman-filter-based estimation of the state of
%                     quadcopter from sensor measurements.
%
% INPUTS
%
%          delt = Measurement update interval, in seconds.
%
%          Tsim = Total time of simulation, in seconds.
%
%      plotFlag = Binary flag determining whether or not the specified
%                 trajectory is plotted.
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
%        xIstar = Nx3 matrix of desired x-axis unit vectors with respect to
%                 the I frame in meters. xIstar(k,:)' is the 3x1 x-axis
%                 unit vector at time tk = tVec(k).
%
%+------------------------------------------------------------------------------+
% Author: Evan Srnka
%+==============================================================================+
% Prepare Minimum Snap pathplanning algorithm
%Polynomial order
W.O = 7;  
%Time to reach each waypoint
Winterval = 2;
W.tVecWp = [0:Winterval:Tsim]; 
%Output sampling period
W.Ts = delt;            

% Conditions on x
% Position x at each waypoint
x_w = [0 1 0 -1 0 1 0 -1 0 1 0 ...
       0 1 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6]; 
% Velocity x at each waypoint
vx_w = zeros(1,length(W.tVecWp));
vx_w(2:end-1) = NaN;
% Acceleration x at each waypoint
ax_w = zeros(1,length(W.tVecWp));
ax_w(2:end-1) = NaN;
% Jerk x at each waypoint
jx_w = zeros(1,length(W.tVecWp));
jx_w(2:end-1) = NaN;   

% Conditions on y
% Position y at each waypoint
y_w = [0 0 1 0 -1 0 1 0 -1 0 1 ...
        0 -1 0 1 0 -1 0 1 0 -1]; 
% Velocity y at each waypoint
vy_w = zeros(1,length(W.tVecWp));
vy_w(2:end-1) = NaN;
% Acceleration y at each waypoint
ay_w = zeros(1,length(W.tVecWp));
ay_w(2:end-1) = NaN;   
% Jerk y at each waypoint
jy_w = zeros(1,length(W.tVecWp));
jy_w(2:end-1) = NaN;  

% Conditions on z
% Position z at each waypoint
z_w = [0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 ...
       3.0 2.0 1.0 2.0 3.0 2.0 1.0 2.0 3.0 2.0];       
% Velocity y at each waypoint
vz_w = zeros(1,length(W.tVecWp));
vz_w(2:end-1) = NaN;   
% Acceleration y at each waypoint
az_w = zeros(1,length(W.tVecWp));
az_w(2:end-1) = NaN;
%Jerk x at each waypoint
jz_w = zeros(1,length(W.tVecWp));
jz_w(2:end-1) = NaN;    

W.rIwp = [x_w' y_w' z_w'];
W.vIwp = [vx_w' vy_w' vz_w'];
W.aIwp = [ax_w' ay_w' az_w'];
W.jIwp = [jx_w' jy_w' jz_w'];

%  Solve minimum snap to populate reference trajectory
R = PathSmoothing(W);
R.xIstar = R.rIstar/norm(R.rIstar);

if(plotFlag)
    % Plot of desired trajectory
    plot3(R.rIstar(:,1), R.rIstar(:,2), R.rIstar(:,3),'Linewidth',2); 
    hold on; grid on;
    axislim = [min(R.rIstar(:,1))-0.5 max(R.rIstar(:,1))+0.5 ...
        min(R.rIstar(:,2))-0.5 max(R.rIstar(:,2))+0.5 ...
        min(R.rIstar(:,3))-0.5 max(R.rIstar(:,3))+0.5];
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
end

end