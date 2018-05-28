% Top-level script for calling simulateQuadrotorControl or
% simulateQuadrotorEstimationAndControl. This version controls the neural
% network version of the simulator

% Use clear all to clear out persistent variables from run to run
clear all; clc;

%Set the random number generator seed for consistency
rng(0);
% Assert this flag to call the full estimation and control simulator;
% otherwise, the control simulator is called
estimationFlag = 1;
% Assert this flag to plot the planned trajectory, otherwise you will only
% see the results
plotFlag = 1;
% Total simulation time, in seconds
Tsim = 40;
% Control update interval, in seconds (corresponds to a 200Hz IMU)
delt = 0.005;

% Construct the time vector, in seconds 
N = floor(Tsim/delt);
tVec=[0:N]'*delt;

% Perform pathplanning
R = pathplan(delt, Tsim, plotFlag);

% Matrix of disturbance forces acting on the body, in N, expressed in I.
S.distMat= 0*3*randn(N+1,3);
% Initial position in m
S.state0.r=[0 0 0]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e=[0 0 0]';
% Initial velocity of body with respect to I, expressed in I, in m/s
S.state0.v=[0 0 0]';
% Initial angular rate of body with respect to I, expressed in B, in rad/s
S.state0.omegaB=[0 0 0]';
% Oversampling factor
S.oversampFact=2;
% Feature locations in the I frame
S.rXIMat = [0,0,1; 0,0,0.5; 1,1,1; 0,0,2; 1,0,1; 1,0,2; 2,0,0; 3,0,0];
% Load quadrotor parameters and constants from files
quadParamsScript;
constantsScript;
sensorParamsScript;
% Load parameters into a structure to be passed into functions
P.quadParams = quadParams; 
P.constants = constants; 
P.sensorParams = sensorParams;

if(estimationFlag)
  [Q, Ms] = simulateQuadrotorEstimationAndControl(R,S,P);
else
    % This version exists for checking state dynamics update. Assumes truth
    % is known at all times.
  Q = simulateQuadrotorControl(R,S,P);
end

% Begin Plotting

% Real-time visualization of quad location
S2.tVec = Q.tVec;
S2.rMat = Q.state.rMat;
S2.eMat = Q.state.eMat;
S2.plotFrequency = 20;
% Set this flag to true to save the visualization as a gif
S2.makeGifFlag = false;
S2.gifFileName = 'neuralnetGif.gif';
% These bounds are hardcoded to match the planned path and may need to
% change depending on path
S2.bounds=5*[-.6 .6 -.6 .6 -0.15 1.1];
visualizeQuad(S2);




