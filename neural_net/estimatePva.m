function [r,v,a] = estimatePva(rTildeMat,delt)
% estimatePva : Estimate the 3-D position, velocity, and acceleration from
%               a time history of position measurements. 
%
% INPUTS
%
% rTildeMat -- Nx3 matrix of position measurements, in meters.
%              rTildeMat(i,:)' is the 3x1 ith position measurement.
%
% delt ------- Measurement sampling time, in seconds
%
%
% OUTPUTS
% 
% r ---------- 3x1 estimate of position that applies at the instant
%              rTildeMat(N,:)' was measured.
%
% v ---------- 3x1 estimate of velocity that applies at the instant
%              rTildeMat(N,:)' was measured.
%
% a ---------- 3x1 estimate of acceleration that applies at the instant
%              rTildeMat(N,:)' was measured.
%
%+------------------------------------------------------------------------------+
% Author:  Todd Humphreys
%+==============================================================================+  

Np = 2;
r = zeros(3,1); v = zeros(3,1); a = zeros(3,1);
[N,~] = size(rTildeMat);
tVec = [0:N-1]'*delt;

if(N < Np + 1)
  error('Polyfit requires at least Np measurement samples.');
end

Px = polyfit(tVec,rTildeMat(:,1),Np);
Py = polyfit(tVec,rTildeMat(:,2),Np);
Pz = polyfit(tVec,rTildeMat(:,3),Np);

% Position estimate 
r(1) = polyval(Px,tVec(end));
r(2) = polyval(Py,tVec(end));
r(3) = polyval(Pz,tVec(end));

% Velocity estimate 
v(1) = polyval(polyder(Px),tVec(end));
v(2) = polyval(polyder(Py),tVec(end));
v(3) = polyval(polyder(Pz),tVec(end));

% Acceleration estimate
a(1) = polyval(polyder(polyder(Px)),tVec(end));
a(2) = polyval(polyder(polyder(Py)),tVec(end));
a(3) = polyval(polyder(polyder(Pz)),tVec(end));


