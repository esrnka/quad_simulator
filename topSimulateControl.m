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
Tsim = 40;
% Control update interval, in seconds
delt = 0.005;
% Time vector, in seconds 
N = floor(Tsim/delt);
tVec=[0:N]'*delt;

% Perform pathplanning
R = pathplan(delt, Tsim, 1);

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
S.rXIMat = [0,0,1; 0,0,0.5; 1,1,1; 0,0,2; 1,0,1; 1,0,2; 2,0,0];
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
%[rXIHat, Re] = estimate3dFeatureLocation(Ms,P);

% Begin Plotting

% Real-time visualization of quad location
S2.tVec = Q.tVec;
S2.rMat = Q.state.rMat;
S2.eMat = Q.state.eMat;
S2.plotFrequency = 20;
S2.makeGifFlag = false;
S2.gifFileName = 'testGif.gif';
S2.bounds=5*[-.6 .6 -.6 .6 -0.15 1.1];
visualizeQuad(S2);

% Plot of norm of lever arm covariance
figure(3); plot(tVec(1:end-1),Q.state.PlB); xlabel('t (s)'); ylabel('Norm of P_{lB}');
text(5, Q.state.PlB(1)/2, ['Final lever arm: ' num2str(Q.state.lB(end,:))]);
% figure(2);clf;
% plot(Q.tVec,Q.state.rMat(:,3)); grid on;
% xlabel('Time (sec)');
% ylabel('Vertical (m)');
% title('Vertical position of CM'); 
% 
% figure(3);clf;
% psiError = unwrap(Q.tVec + pi - Q.state.eMat(:,3));
% meanPsiErrorInt = round(mean(psiError)/2/pi);
% plot(Q.tVec,psiError - meanPsiErrorInt*2*pi);
% grid on;
% xlabel('Time (sec)');
% ylabel('\Delta \psi (rad)');
% title('Yaw angle error');
% 
% figure(5);clf;
% plot(Q.state.rMat(:,1), Q.state.rMat(:,2)); 
% axis equal; grid on;
% xlabel('X (m)');
% ylabel('Y (m)');
% title('Horizontal position of CM');



