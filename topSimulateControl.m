% Top-level script for calling simulateQuadrotorControl or
% simulateQuadrotorEstimationAndControl

% Need to use clear all to clear out persistent variables from run to run
clear all; clc;

%Set the random number generator seed
rng(0);
% Assert this flag to call the full estimation and control simulator;
% otherwise, the control simulator is called
estimationFlag = 1;
% Total simulation time, in seconds
Tsim = 10;
% Control update interval, in seconds
delt = 0.005;
% Time vector, in seconds 
N = floor(Tsim/delt);
tVec=[0:N-1]'*delt;

% Prepare Minimum Snap pathplanning algorithm
W.O = 7;                %Polynomial order
W.tVecWp = [0 2 2 2 2];   %Time to reach each waypoint
W.Ts = delt;            %Output sampling period

%Conditions on x
x_w = [0 1 2 3 4];   %Position x at each waypoint
vx_w = [0 NaN NaN NaN 0];    %Velocity x at each waypoint
ax_w = [0 NaN NaN NaN 0];    %Acceleration x at each waypoint
jx_w = [0 NaN NaN NaN 0];    %Jerk x at each waypoint

%Conditions on y
y_w = [0 0 1 2 1];       %Position y at each waypoint
vy_w = [0 NaN NaN NaN 0];    %Velocity y at each waypoint
ay_w = [0 NaN NaN NaN 0];    %Acceleration y at each waypoint
jy_w = [0 NaN NaN NaN 0];    %Jerk y at each waypoint

%Conditions on z
z_w = [0 1 1 1 1];       %Position z at each waypoint
vz_w = [0 NaN NaN NaN 0];    %Velocity y at each waypoint
az_w = [0 NaN NaN NaN 0];    %Acceleration y at each waypoint
jz_w = [0 NaN NaN NaN 0];    %Jerk x at each waypoint

W.rIwp = [x_w' y_w' z_w'];
W.vIwp = [vx_w' vy_w' vz_w'];
W.aIwp = [ax_w' ay_w' az_w'];
W.jIwp = [jx_w' jy_w' jz_w'];

%  Solve minimum snap to populate reference trajectory
R = PathSmoothing(W);
R.xIstar = -R.rIstar/norm(R.rIstar);

disp("Path Planning Complete")

% Matrix of disturbance forces acting on the body, in N, expressed in I.
S.distMat= 0*3*randn(N-1,3);
% Initial position in m
S.state0.r=[x_w(1) y_w(1) z_w(1)]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e=[0 0 0]';
% Initial velocity of body with respect to I, expressed in I, in m/s
S.state0.v=[vx_w(1) vy_w(1) vz_w(1)]';
% Initial angular rate of body with respect to I, expressed in B, in rad/s
S.state0.omegaB=[0 0 0]';
% Oversampling factor
S.oversampFact=2;
% Feature locations in the I frame
S.rXIMat = [0,0,1; 0,0,0.5; 1,1,1];
% Quadrotor parameters and constants
quadParamsScript;
constantsScript;
sensorParamsScript;
P.quadParams = quadParams; 
P.constants = constants; 
P.sensorParams = sensorParams;

if(estimationFlag)
  [Q, Ms] = simulateQuadrotorEstimationAndControl(R,S,P);
else
  Q = simulateQuadrotorControl(R,S,P);
end

% Estimate the feature location
[rXIHat, Re] = estimate3dFeatureLocation(Ms,P);

S2.tVec = Q.tVec;
S2.rMat = Q.state.rMat;
S2.eMat = Q.state.eMat;
S2.plotFrequency = 20;
S2.makeGifFlag = false;
S2.gifFileName = 'testGif.gif';
S2.bounds=2.5*[-1 1 -1 1 -0.1 1];
visualizeQuad(S2);

figure(2);clf;
plot(Q.tVec,Q.state.rMat(:,3)); grid on;
xlabel('Time (sec)');
ylabel('Vertical (m)');
title('Vertical position of CM'); 

figure(3);clf;
psiError = unwrap(n*Q.tVec + pi - Q.state.eMat(:,3));
meanPsiErrorInt = round(mean(psiError)/2/pi);
plot(Q.tVec,psiError - meanPsiErrorInt*2*pi);
grid on;
xlabel('Time (sec)');
ylabel('\Delta \psi (rad)');
title('Yaw angle error');

figure(5);clf;
plot(Q.state.rMat(:,1), Q.state.rMat(:,2)); 
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Horizontal position of CM');



