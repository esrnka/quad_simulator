function [zk] = h_meas(xk,wk,RBIBark,vIMat,P)
% h_meas : Measurement model for quadcopter.
%
% INPUTS
%
% xk --------- 15x1 state vector at time tk, defined as 
% 
%              xk = [rI', vI', e', ba', bg']'
%
%              where all corresponding quantities are identical to those
%              defined for E.statek in stateEstimatorUKF.m and where e is the
%              3x1 error Euler angle vector defined such that for an estimate
%              RBIHat of the attitude, the true attitude is RBI =
%              dRBI(e)*RBIHat, where dRBI(e) is the DCM formed from the error
%              Euler angle vector e.
%
% wk --------- nz-by-1 measurement noise vector at time tk, defined as
%
%              wk = [wPIk', wCIk', w1B', w2B', ..., wNfB']'
%
%              where nz = 6 + Nf*3, and where all 3x1 noise vectors represent
%              additive noise on the corresponding measurements.
%
% RBIBark ---- 3x3 attitude matrix estimate at time tk.
%
% vIMat ------ Nf-by-3 matrix of coordinates of visual features in the
%              simulation environment, expressed in meters in the I
%              frame. vIMat(i,:)' is the 3x1 vector of coordinates of the ith
%              feature.
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
% zk --------- nz-by-1 measurement vector at time tk, defined as
%
%              zk = [rPItilde', rCItildeu', v1Btildeu', ..., vNfBtildeu']'
%
%              where rPItilde is the 3x1 measured position of the primary
%              antenna in the I frame, rCItildeu is the 3x1 measured unit
%              vector pointing from the primary to the secondary antenna,
%              expressed in the I frame, and viBtildeu is the 3x1 unit vector,
%              expressed in the body frame, pointing toward the ith 3D
%              feature.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Todd Humphreys
%+==============================================================================+  


rIk = xk(1:3);
ek = xk(7:9);
rPB = P.sensorParams.rA(:,1);
rCB = P.sensorParams.rA(:,2) - P.sensorParams.rA(:,1);
rCBu = rCB/norm(rCB);
RBIk = euler2dcm(ek)*RBIBark;
hk = [rIk + RBIk'*rPB; RBIk'*rCBu];
[Nf,~] = size(vIMat);
for ii=1:Nf
  hk = [hk; RBIk*vIMat(ii,:)'];
end

zk = hk + wk;
