function [rXIHat,Re] = estimate3dFeatureLocation(M,P)
% estimate3dFeatureLocation : Estimate the 3D coordinates of a feature
%                             point seen by two or more cameras with known
%                             pose.    
%
%
% INPUTS
%
% M ---------- Structure with the following elements:
%
%       rxArray = Nx1 cell array of measured positions of the feature point
%                 projection on the camera's image plane, in pixels.
%                 rxArray{i} is the 2x1 vector of coordinates of the feature
%                 point as measured by the ith camera.  To ensure the
%                 estimation problem is observable, N must satisfy N >= 2 and
%                 at least two cameras must be non-colinear.
%
%      RCIArray = Nx1 cell array of camera-to-I-frame attitude matrices.
%                 RCIArray{i} is the 3x3 attitude matrix corresponding to the
%                 measurement rxArray{i}.
%
%       rcArray = Nx1 cell array of camera center positions.  rcArray{i} is
%                 the 3x1 position of the camera center corresponding to the
%                 measurement rxArray{i}, expressed in the I frame in meters.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
% OUTPUTS
% rXIHat -------- 3x1 estimated location of the feature point expressed in I
%                 in meters.
%
% Re ------------ 3x3 error covariance matrix for the estimate rxIHat.
%
%+------------------------------------------------------------------------------+
% Author:  Evan Srnka
%+==============================================================================+  

%Determine number of features 
N = size(M.rxArray,1);
%Create H and R matrices
H = [];
R = [];
%Create intrinsic paramter matrix k 
k = P.sensorParams.K; 
%Extract the measurments
for i = 1:N
    Proj = k * [M.RCIArray{i},-M.RCIArray{i} * M.rcArray{i}];
    X = [P.sensorParams.pixelSize * M.rxArray{i}];
    %Projection Row vectors
    P1 = Proj(1,:);
    P2 = Proj(2,:);
    P3 = Proj(3,:);    
    %Stack the H matrix with the new feature
    H = [H; 
         X(1)*P3 - P1; 
         X(2)*P3 - P2 ];    
    % Stack the R matrix with the new noise
    R = blkdiag(R, P.sensorParams.pixelSize^2 * P.sensorParams.Rc);
end
%Find Hr and z
Hr = H(:,1:3);
z = -1 * H(:,end);
% Least Squares Method
rXIHat = inv(Hr'*inv(R) * Hr)*Hr'*inv(R)*z;
Re = inv(Hr'*inv(R)*Hr);
