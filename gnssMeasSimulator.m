function [rPtilde,rStilde,rCtilde] = gnssMeasSimulator(S,P)
% gnssMeasSimulator : Simulates GNSS measurements for quad.
%
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
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude of B frame wrt I frame
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
%
% OUTPUTS
%
% rPtilde----- 3x1 GNSS-measured position of the quad's primary GNSS antenna,
%              in ECEF coordinates relative to the reference antenna, in
%              meters.
%
% rStilde ---- 3x1 GNSS-measured position of the quad's secondary GNSS
%              antenna, in ECEF coordinates relative to the reference antenna,
%              in meters.
%
% rCtilde ---- 3x1 GNSS-measured position of secondary GNSS antenna, in ECEF
%              coordinates relative to the primary antenna, in meters.
%              rCtilde is constrained to satisfy norm(rCtilde) = b, where b is
%              the known distance between the two antennas.
%
%+------------------------------------------------------------------------------+
% Author: Evan Srnka
%+==============================================================================+  

RIG = Recef2enu(P.sensorParams.r0G);
RPa = chol(RIG'*P.sensorParams.RpL*RIG);
RIB = S.statek.RBI'; 
rPI = S.statek.rI + RIB*P.sensorParams.rA(:,1);
rPG = RIG'*rPI + P.sensorParams.r0G;
rP = rPG - P.sensorParams.rRefG;
rPtilde = rP + RPa'*randn(3,1);
rSI = S.statek.rI + RIB*P.sensorParams.rA(:,2);
rSG = RIG'*rSI + P.sensorParams.r0G;
rS = rSG - P.sensorParams.rRefG;
rStilde = rS + RPa'*randn(3,1);
rC = rSG - rPG; 
rCu = rC/norm(rC);
% Add an epsilon along the diagonal to ensure RC is positive definite so we
% can apply Cholesky decomposition
epsilon = 1e-8;
RC = P.sensorParams.sigmaC^2*(eye(3)-rCu*rCu') + epsilon*eye(3);
wC = chol(RC)'*randn(3,1);
rCtilde = rC + wC;
rCtilde = norm(rC)*rCtilde/norm(rCtilde);

