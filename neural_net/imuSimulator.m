function [fB,omegaBtilde] = imuSimulator(S,P)
% imuSimulator : Simulates IMU measurements of specific force and
%                body-referenced angular rate.
%
% INPUTS
%
% S ---------- Structure with the following elements:
%
%        statek = State of the quad at tk, expressed as a structure with the
%                 following elements:
%                   
%                  rI = 3x1 position of CM in the I frame, in meters
%
%                  vI = 3x1 velocity of CM with respect to the I frame and
%                       expressed in the I frame, in meters per second.
%
%                  aI = 3x1 acceleration of CM with respect to the I frame and
%                       expressed in the I frame, in meters per second^2.
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude of B frame wrt I frame
%
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
%           omegaBdot = 3x1 time derivative of omegaB, in radians per
%                       second^2.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
% OUTPUTS
%
% fB --------- 3x1 specific force measured by the IMU's 3-axis accelerometer
%
% omegaBtilde  3x1 angular rate measured by the IMU's 3-axis rate gyro
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+
persistent ba bg

if(isempty(ba))
  % Set ba's initial value
  QbaSteadyState = P.sensorParams.Qa2/(1 - P.sensorParams.alphaa^2);
  ba = (chol(QbaSteadyState))'*randn(3,1);
end
if(isempty(bg))
  % Set bg's initial value
  QbgSteadyState = P.sensorParams.Qg2/(1 - P.sensorParams.alphag^2);
  bg = (chol(QbgSteadyState))'*randn(3,1);
end

% Propagate biases
ba = P.sensorParams.alphaa*ba + chol(P.sensorParams.Qa2)'*randn(3,1);
bg = P.sensorParams.alphag*bg + chol(P.sensorParams.Qg2)'*randn(3,1);

% Define misalignment and scale factor error matrix (from accelgen.m)
M = [1.0682, 0.0014, 0.0792; -0.0208, 0.9950, 0.0436; -0.0921, ...
    0.0232, 1.0936];

% IMU measurements
fB = M*S.statek.RBI*(S.statek.aI + [0;0;P.constants.g]) + ...
     cross(S.statek.omegaB,cross(S.statek.omegaB,P.sensorParams.lB)) + ...
     cross(S.statek.omegaBdot,P.sensorParams.lB) + ba + ... 
     chol(P.sensorParams.Qa)'*randn(3,1);

omegaBtilde = M*S.statek.omegaB + bg + chol(P.sensorParams.Qg)'*randn(3,1);

