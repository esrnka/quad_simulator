function [rx] = hdCameraSimulator(rXI,S,P)
% hdCameraSimulator : Simulates feature location measurements from the
%                     quad's high-definition camera. 
%
%
% INPUTS
%
% rXI -------- 3x1 location of a feature point expressed in I in meters.
%
% S ---------- Structure with the following elements:
%
%        statek  = State of the quad at tk, expressed as a structure with the
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
% OUTPUTS
%
% rx --------- 2x1 measured position of the feature point projection on the
%              camera's image plane, in pixels.  If the feature point is
%              not visible to the camera (the ray from the feature to the
%              camera center never intersects the image plane), then rx is
%              an empty matrix.
%
%+------------------------------------------------------------------------------+
% Author:  Evan Srnka
%+==============================================================================+  

rx = [];

% Construct transformation matrix P = K*[RCI,t]
RCI = P.sensorParams.RCB * S.statek.RBI;
t = -RCI*(S.statek.rI + S.statek.RBI'*P.sensorParams.rc);
Pk = P.sensorParams.K*[RCI t];

% Generate projection
x = Pk*[rXI;1];
xc = x(1:2)/x(3); 

% Check if point is within image plane
if(x(3) > 0 && abs(xc(1)) < P.sensorParams.imagePlaneSize(1)/2 && ...
   abs(xc(2)) < P.sensorParams.imagePlaneSize(2)/2)
  
  Rac = chol(P.sensorParams.Rc);
  rx = (1/P.sensorParams.pixelSize)*xc + Rac'*randn(2,1);
end
