function [eak] = voltageConverter(Fk,NBk,P)
% trajectoryController : Controls quadcopter toward a reference trajectory.
%
%
% INPUTS
%
% Fk --------- Commanded total thrust at time tk, in Newtons.
%
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%
% OUTPUTS
%
% eak -------- Commanded 4x1 voltage vector to be applied at time tk, in
%              volts. eak(i) is the voltage for the ith motor.
%
%+------------------------------------------------------------------------------+
% Author:  Evan Srnka
%+==============================================================================+  

% Determine maximum force that can be applied by any rotor
omegaMax = min(P.quadParams.cm) * P.quadParams.eamax; 
FMax = min(P.quadParams.kF)*(omegaMax^2);

% Populate conversion matrix
kTVec = P.quadParams.kN./P.quadParams.kF;
G = [ones(1,4); 
     P.quadParams.rotor_loc(2,:);
     -P.quadParams.rotor_loc(1,:);
     -(kTVec').*P.quadParams.omegaRdir];
GInv = inv(G);

% Apply non-negative force and saturation constraints
beta = 0.95;
Fks = max(min(Fk,4*beta*FMax), 4*(1-beta)*FMax);
alpha = 1;
FVec = GInv*[Fks;alpha*NBk];
while(max(FVec) > FMax || min(FVec) < 0)
  alpha = 0.99*alpha;
  FVec = GInv*[Fks;alpha*NBk];
end

% Convert force to voltage
omegaVec = sqrt(FVec./P.quadParams.kF);
eak = omegaVec./P.quadParams.cm;


