function [E] = stateEstimator(S,M,P)
% stateEstimator : Estimates state of quadcopter from sensor measurements.  It
%                  is assumed that during the first 1 second of sensor
%                  measurements the quad is at rest in the I frame with with
%                  [0;0;1] = RBI*[0;0;1].
%
%
% INPUTS
%
% S ---------- Structure with the following elements:
%
%        rXIMat = Nf-by-3 matrix of coordinates of visual features in the
%                 simulation environment, expressed in meters in the I
%                 frame. rXIMat(i,:)' is the 3x1 vector of coordinates of
%                 the ith feature.  
%
%          delt = Measurement update interval, in seconds.
%
% M ---------- Structure with the following elements:
%
%            tk = Time at which all measurements apply, in seconds.
%
%       rPtilde = 3x1 GNSS-measured position of the quad's primary GNSS
%                 antenna, in ECEF coordinates relative to the reference
%                 antenna, in meters.
%
%       rStilde = 3x1 GNSS-measured position of the quad's secondary GNSS
%                 antenna, in ECEF coordinates relative to the reference
%                 antenna, in meters.
%
%       rCtilde = 3x1 GNSS-measured position of secondary GNSS antenna, in
%                 ECEF coordinates relative to the primary antenna, in meters.
%                 rC is constrained to satisfy norm(rC) = b, where b is the
%                 known distance between the two antennas.
%
%         rxMat = Nf-by-2 matrix of measured positions of feature point
%                 projections on the camera's image plane, in
%                 pixels. rxMat(i,:)' is the 2x1 image position measurement of
%                 the 3D feature in rXIMat(i,:)'.  If the ith feature point is
%                 not visible to the camera (the ray from the feature to the
%                 camera center never intersects the image plane), then
%                 rxMat(i,:) = [NaN, NaN].               
%
%            fB = 3x1 specific force measured by the IMU's 3-axis
%                 accelerometer
%
%   omegaBtilde = 3x1 angular rate measured by the IMU's 3-axis rate gyro
%
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%  sensorParams = Structure containing sensor parameters, as defined in
%                 sensorParamsScript.m
%
%
% OUTPUTS
%
% E ---------- Structure with the following elements:
%
%        statek = Estimated state of the quad at tk, expressed as a structure
%                 with the following elements.  During the initial 1-second
%                 calibration period, statek is an empty matrix.
%                   
%                  rI = 3x1 position of CM in the I frame, in meters
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude
%
%                  vI = 3x1 velocity of CM with respect to the I frame and
%                       expressed in the I frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
%                  ba = 3x1 bias of accelerometer expressed in the
%                       accelerometer frame in meter/second^2.
%
%                  bg = 3x1 bias of rate gyro expressed in the body frame in
%                       rad/sec.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+  

%----- Setup
% Averaging interval, in seconds
tAvg = 0.02;
NAvg = floor(tAvg/S.delt);
% Polyfit interval, in seconds
tPoly = 0.05;
NPoly = floor(tPoly/S.delt);
% P-to-S vector in body frame
rCB = P.sensorParams.rA(:,2) - P.sensorParams.rA(:,1);
rCBu = rCB/norm(rCB);
e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1];
RIG = Recef2enu(P.sensorParams.r0G);
% Length of initial calibration interval, in seconds
tCal = 1;

persistent statek D calibrated 
if(isempty(statek))
  % Set initial state estimate values
  statek.rI = zeros(3,1);
  statek.RBI = eye(3,3);
  statek.vI = zeros(3,1);
  statek.omegaB = zeros(3,1);
  statek.ba = zeros(3,1);
  statek.bg = zeros(3,1); 
  D.tVec = [];
  D.rPItildeMat = [];
  D.rSItildeMat = [];
  D.rCItildeMat = [];
  D.fBMat = [];
  D.omegaBtildeMat = [];
  calibrated = false;
end

%----- Store incoming measurements
D.tVec = [D.tVec; M.tk];
% Convert GNSS measurements to I frame
rPItilde = RIG*(M.rPtilde + P.sensorParams.rRefG - P.sensorParams.r0G);
rSItilde = RIG*(M.rStilde + P.sensorParams.rRefG - P.sensorParams.r0G);
rCItilde = RIG*M.rCtilde;
D.rPItildeMat = [D.rPItildeMat; rPItilde'];
D.rSItildeMat = [D.rSItildeMat; rSItilde'];
D.rCItildeMat = [D.rCItildeMat; rCItilde'];
D.fBMat = [D.fBMat; M.fB'];
D.omegaBtildeMat = [D.omegaBtildeMat; M.omegaBtilde'];

%----- Calibrate
% The quad is guaranteed to be at rest in the I frame with ZI = RBI'*zB for
% the first 1 second of measurements.  We exploit this interval to obtain
% initial estimates of bg, ba, and RBI.
if(M.tk - D.tVec(1) >= tCal && ~calibrated)
  statek.ba = mean(D.fBMat,1)' - [0;0;P.constants.g];
  statek.bg = mean(D.omegaBtildeMat,1)';
  % Get initial estimate of RBI
  rCIHat = mean(D.rCItildeMat,1)';
  rCIHatu = rCIHat/norm(rCIHat);
  vIMat = [rCIHatu'; e3'];
  vBMat = [rCBu'; e3'];
  aVec = ones(2,1);
  statek.RBI = wahbaSolver(aVec,vIMat,vBMat);
  calibrated = true;
end
if(~calibrated)
  E.statek = [];
  return;
end

%----- Estimate angular rate and acceleration
NPolyOmega = NPoly;
[omegaBtildef,omegaBtildeDotf,~] = ...
    estimatePva(D.omegaBtildeMat(end - NPolyOmega:end,:),S.delt);
statek.omegaB = omegaBtildef - statek.bg;
omegaBdotHat = omegaBtildeDotf;

%----- Estimate CM position and velocity
[rPItildef,vPItildef,aPItildef] = ...
    estimatePva(D.rPItildeMat(end - NPoly:end,:),S.delt);
rPB = P.sensorParams.rA(:,1);
statek.rI = rPItildef - (statek.RBI')*rPB;
statek.vI = vPItildef - (statek.RBI')*cross(statek.omegaB,rPB);
aIHat = aPItildef - ...
        (statek.RBI')*cross(statek.omegaB,cross(statek.omegaB,rPB)) - ...
        (statek.RBI')*cross(omegaBdotHat,rPB);

%----- Estimate RBI
% Estimate rC vector in I frame
[rCItildef,~,~] = estimatePva(D.rCItildeMat(end - NPoly:end,:),S.delt);
rCIHatu = rCItildef/norm(rCItildef);
% Estimate anti-gravity vector in B frame
fBf = mean(D.fBMat(end - NAvg:end,:),1)';
lB = P.sensorParams.lB;
gBHat = fBf - (statek.RBI)*aIHat - cross(omegaBdotHat,lB) - ...
        cross(statek.omegaB,cross(statek.omegaB,lB)) - statek.ba;
gBHatu = gBHat/norm(gBHat);
gIu = [0;0;1];
% Obtain pairs of vectors pointing to 3D feature points in B and I frames
vIMat = []; vBMat = []; aVec = [];
[Nft,~] = size(M.rxMat);
for ii=1:Nft
  if(~isnan(M.rxMat(ii,1)))
    viC = [P.sensorParams.pixelSize*M.rxMat(ii,:)';P.sensorParams.f];
    viCu = viC/norm(viC);
    viBu = (P.sensorParams.RCB')*viCu;
    % Position of camera center in I frame
    rcIHat = statek.rI + (statek.RBI')*P.sensorParams.rc;
    viI = S.rXIMat(ii,:)' - rcIHat;
    viIu = viI/norm(viI);
    vIMat = [vIMat; viIu'];
    vBMat = [vBMat; viBu'];
    aVec = [aVec; 1];
  end
end
% Add rC and anti-gravity vector to vectors used for attitude estimation
vIMat = [vIMat; rCIHatu'; gIu'];
vBMat = [vBMat; rCBu'; gBHatu'];
% Weight the rC vector as much as the camera-derived vectors, but the gravity
% vector less, since it's based on a noisy IMU and other estimates
aVec = [aVec; 1; 0.01];
% Propagate the prior estimate of RBI forward using the estimated angular rate
% and basic Euler integration to get a new estimate of RBI at the current time
RBInew = statek.RBI - S.delt*crossProductEquivalent(statek.omegaB)*statek.RBI;
% Extract three pseudo-measurement vectors from RBInew
e1B = RBInew*e1; e2B = RBInew*e2; e3B = RBInew*e3;
% Combine the pseudo-measurement vectors with the other measured vectors to be
% used for attitude estimation and weight the pseudo-measurement vectors
% appropriately
vIMat = [vIMat; e1'; e2'; e3'];
vBMat = [vBMat; e1B'; e2B'; e3B'];
aVec = [aVec; 0.5; 0.5; 0.5];
[statek.RBI] = wahbaSolver(aVec,vIMat,vBMat);

if(0)
% Testing section
gBHatuTrue = statekTrue.RBI*gIu;
(180/pi)*acos(gBHatuTrue'*gBHatu)
%[rPITrue, rPItildef]
[statekTrue.aI, aIHat]
%[statekTrue.omegaB, statek.omegaB]
[statekTrue.omegaBdot, omegaBdotHat]
dcm2euler(statekTrue.RBI' * statek.RBI)*180/pi
pause;
clc;
end

E.statek = statek;