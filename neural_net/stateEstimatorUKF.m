function [E] = stateEstimatorUKF(S,M,P)
% stateEstimatorUKF : Unscented-Kalman-filter-based estimation of the state of
%                     quadcopter from sensor measurements.  
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
%           lB = lever arm
%
% OUTPUTS
%
% E ---------- Structure with the following elements:
%
%        statek = Estimated state of the quad at tk, expressed as a structure
%                 with the following elements. 
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
%            Pk = 15x15 error covariance matrix for the estimator state
%                 vector xk that applies at time tk, which is defined as 
% 
%                 xk = [rI', vI', e', ba', bg']'
%
%                 where all corresponding quantities are identical to those
%                 defined for E.statek and where e is the 3x1 error Euler
%                 angle vector defined such that for an estimate RBIHat of the
%                 attitude, the true attitude is RBI = dRBI(e)*RBIHat, where
%                 dRBI(e) is the DCM formed from the error Euler angle vector
%                 e, expressed in radians.
%
%+------------------------------------------------------------------------------+
% Author: Todd Humphreys, Evan Srnka
%+==============================================================================+  

%----- Setup
persistent xBark PBark RBIBark
RIG = Recef2enu(P.sensorParams.r0G);                                        %Rotation matrix from G to I
rCB = P.sensorParams.rA(:,2) - P.sensorParams.rA(:,1);                      %Body frame vector from secondary antenna to primary
rCBu = rCB/norm(rCB);                                                       %Unit normalized version of rCB
rPB = P.sensorParams.rA(:,1);                                               %Location of primary antenna in body frame
e3 = [0;0;1];                                                               %z axis vector
epsilon = 1e-8;                                                             %
nx = 15;                                                                    % size of state vector
nv = 12;                                                                    % size of process noise vector
alphaUKF = 1e-3;                                                            % UKF tuning parameter
betaUKF = 2;                                                                % UKF tuning parameter
kappaUKF = 0;                                                               % UKF tuning parameter
lambda_p = alphaUKF^2*(kappaUKF + nx + nv) - nx - nv;                       % UKF transform parameter
c_p = sqrt(nx+nv+lambda_p);                                                 % UKF transform parameter
Wm0_p = lambda_p/(nx + nv + lambda_p);                                      % 0th mean weight
Wmi_p = 1/(2*(nx + nv + lambda_p));                                         % ith mean weight
Wc0_p = Wm0_p + 1 - alphaUKF^2 + betaUKF;                                   % 0th covariance weight
Wci_p = Wmi_p;                                                              % ith covariance weight

%----- Convert GNSS measurements to I frame
rPItilde = RIG*(M.rPtilde + P.sensorParams.rRefG - P.sensorParams.r0G);     %Location of primary antenna in inertial frame
rSItilde = RIG*(M.rStilde + P.sensorParams.rRefG - P.sensorParams.r0G);     %Location of secondary antenna in inertial frame
rCItilde = RIG*M.rCtilde;                                                   %Inertial frame vector between primary and secondary antenna
rCItildeu = rCItilde/norm(rCItilde);                                        %Normalized version of rCItilde

%----- Initialize state estimate on first call to this function
if(isempty(xBark))
  vIMat = [rCItildeu'; e3'];                                                % Pair of inertial frame vectors
  vBMat = [rCBu'; e3'];                                                     % Pair of body frame vectors
  aVec = ones(2,1);                                                         % Weight vector 
  RBIBark = wahbaSolver(aVec,vIMat,vBMat);                                  % Initial estimate of RBI matrix
  rIBark = rPItilde - RBIBark'*rPB;                                         % Inertial frame position
  xBark = [rIBark; zeros(12,1)];                                            % Initial state vector estimate
  QbaSteadyState = P.sensorParams.Qa2/(1 - P.sensorParams.alphaa^2);        % Steady state error covariance matrix of accelerometer bias
  QbgSteadyState = P.sensorParams.Qg2/(1 - P.sensorParams.alphag^2);        % Steady state error covariance matrix of gyro bias
  PBark = diag([0.0025*ones(3,1); 0.001*ones(3,1); 0.01*ones(3,1); ...      % Initial estimate of state covariance matrix
                diag(QbaSteadyState); diag(QbgSteadyState)]);  
end

%----- Assemble measurements
zk = [rPItilde; rCItildeu];                                                 % Form measurement vector from I frame GPS
% Obtain pairs of vectors pointing to 3D feature points in B and I frames
vIMat = []; vBMat = []; RcBCellArray = {}; jj = 1; Nf = 0;
[Nft,~] = size(M.rxMat);                                                    % Number of features = size of camera measurement matrix
for ii=1:Nft
  if(~isnan(M.rxMat(ii,1)))                                                 % If feature is seen by camera
    viCtilde = [P.sensorParams.pixelSize*M.rxMat(ii,:)';P.sensorParams.f];  % Form camera frame measurement estimate
    norm_viCtilde = norm(viCtilde);                                         % Norm this vector
    viCtildeu = viCtilde/norm_viCtilde;                                     % Convert to a unity vector
    viBtildeu = (P.sensorParams.RCB')*viCtildeu;                            % Find the body frame unit vector
    sigmac = sqrt(P.sensorParams.Rc(1,1))*...
        P.sensorParams.pixelSize/norm_viCtilde;                             % Find camera variance
    RcC = sigmac^2*(eye(3) - viCtildeu*viCtildeu') + epsilon*eye(3);        % Form inertial camera covariance matrix
    % Error covariance matrix for the unit vector viBtildeu
    RcB = P.sensorParams.RCB'*RcC*P.sensorParams.RCB;                       % Form body camera covariance matrix
    RcBCellArray{jj} = RcB;                                                 % Store body camera covariance in a cell
    % Position of camera center in I frame
    rcI = xBark(1:3) + RBIBark'*P.sensorParams.rc;                          % Find estiamte of camera center in I frame
    viI = S.rXIMat(ii,:)' - rcI;                                            % Find vector between feature and camera center
    viIu = viI/norm(viI);                                                   % Norm this vector
    vIMat = [vIMat; viIu'];                                                 % Add norm of this vector to measurement matrix
    zk = [zk; viBtildeu];                                                   % Add body frame unit vector to measurement matrix
    jj = jj + 1; Nf = Nf + 1;    
  end
end

%----- Perform measurement update
% Form measurement error covariance matrix Rk
RPI = P.sensorParams.RpL;                                                   % GPS error covariance matrix in ENU
rCIu = RBIBark'*rCBu;                                                       % Estimated constraint vector in inertial
RCI = P.sensorParams.sigmaC^2*(eye(3)-rCIu*rCIu') + epsilon*eye(3);         % Covariance matrix for GPS constraint
Rk = blkdiag(RPI,RCI);                                                      % Measurement error covariance matrix
for ii=1:Nf
  Rk = blkdiag(Rk,RcBCellArray{ii});                                        % Append the feature covariance matrices
end
nz = length(zk);                                                            % Length of measurement vector
lambda_u = alphaUKF^2*(kappaUKF + nx + nz) - nx - nz;                       % UKF transform term
c_u = sqrt(nx+nz+lambda_u);                                                 % UKF transform term
Wm0_u = lambda_u/(nx + nz + lambda_u);                                      % 0th mean weight
Wmi_u = 1/(2*(nx + nz + lambda_u));                                         % ith mean weight
Wc0_u = Wm0_u + 1 - alphaUKF^2 + betaUKF;                                   % 0th covariance weight
Wci_u = Wmi_u;                                                              % ith covariance weight
% Form augmented a priori state and error covariance matrix
xBarAugk = [xBark; zeros(nz,1)];                                            % Augment the predicted state with measurement noise
PBarAugk = blkdiag(PBark,Rk);                                               % Augment the predicted covariance with measurement covariance
SxBar = chol(PBarAugk)';                                                    % Find square root of augmented covariance matrix
% Assemble sigma points and push these through the measurement function
sp0 = xBarAugk;                                                             % Create sigma points equal to augmented state        
spMat = zeros(nx+nz, 2*(nx+nz));                                            % Create zero matrix for sigma points
zpMat = zeros(nz,2*(nx+nz));                                                % Create zero matrix for measurements
zp0 = h_meas(sp0(1:nx),sp0(nx+1:end),RBIBark,vIMat,P);                      % Call the measurement function
for ii=1:2*(nx+nz)
  jj = ii; pm = 1;
  if(ii > (nx + nz)) jj = ii - nx - nz; pm = -1; end
  spMat(:,ii) = sp0 + pm*c_u*SxBar(:,jj);                                   % Sample positive and negative sigma points
  zpMat(:,ii) = h_meas(spMat(1:nx,ii),spMat(nx+1:end,ii),RBIBark,vIMat,P);  % Push sigma points through the measurement function
end
% Recombine sigma points
zBark = sum([Wm0_u*zp0, Wmi_u*zpMat],2);                                    % Calculate predicted measurement from sigma points
Pzz = Wc0_u*(zp0 - zBark)*(zp0 - zBark)';                                   % Calculate predicted measurement covariance
Pxz = Wc0_u*(sp0(1:nx) - xBark)*(zp0 - zBark)';                             % Calculate predicted measurement cross covariance with state
for ii=1:2*(nx+nz)
  Pzz = Pzz + Wci_u*(zpMat(:,ii) - zBark)*(zpMat(:,ii) - zBark)';           % Calculate predicted measurment covariance
  Pxz = Pxz + Wci_u*(spMat(1:nx,ii) - xBark)*(zpMat(:,ii) - zBark)';        % Calculate predicted measurement cross covariacne with state
end
% Perform LMMSE measurement update
PzzInv = inv(Pzz); 
xHatk = xBark + Pxz*PzzInv*(zk - zBark);                                    % LMMSE measurement state update
Pk = PBark - Pxz*PzzInv*Pxz';                                               % LMMSe measurement covariance update

%----- Package output state
statek.rI = xHatk(1:3); statek.vI = xHatk(4:6); ek = xHatk(7:9);
statek.RBI = euler2dcm(ek)*RBIBark;
statek.ba = xHatk(10:12); statek.bg = xHatk(13:15);
statek.omegaB = M.omegaBtilde - statek.bg;
PkDiag = diag(Pk);
E.statek = statek;
E.Pk = Pk;
 
if(0)
% Testing section
tk = M.tk
[statek.rI - statekTrue.rI]
[statek.vI - statekTrue.vI]
dcm2euler(statekTrue.RBI' * statek.RBI)*180/pi
sqrt(PkDiag(7:9))*180/pi
pause;
clc;
end

%----- Propagate state to time tkp1
RBIHatk = statek.RBI;                                                       % Set the best estimate of the RBI matrix
xHatk(7:9) = zeros(3,1);                                                    % Estimate of euler angles to zero
Qk = blkdiag(P.sensorParams.Qg,P.sensorParams.Qg2,...
  P.sensorParams.Qa,P.sensorParams.Qa2);                                    % Covariance matrix of IMU
xHatAugk = [xHatk; zeros(nv,1)];                                            % Augment state vector with process noise
PAugk = blkdiag(Pk,Qk);                                                     % Augment covariance vector with imu covariance
Sx = chol(PAugk)';                                                          % Cholesky square root of augmented covariance
% Assemble sigma points and push these through the dynamics function
sp0 = xHatAugk;
xpMat = zeros(nx,2*(nx+nv));
uk = [M.omegaBtilde;M.fB];
xp0 = f_dynamics(sp0(1:nx),uk,sp0(nx+1:end),S.delt,RBIHatk,P);
for ii=1:2*(nx+nv)
  jj = ii; pm = 1;
  if(ii > (nx + nv)) jj = ii - nx - nv; pm = -1; end
  spii = sp0 + pm*c_p*Sx(:,jj);
  xpMat(:,ii) = f_dynamics(spii(1:nx),uk,spii(nx+1:end),S.delt,RBIHatk,P);
end
% Recombine sigma points
xBarkp1 = sum([Wm0_p*xp0, Wmi_p*xpMat],2);
PBarkp1 = Wc0_p*(xp0 - xBarkp1)*(xp0 - xBarkp1)';
for ii=1:2*(nx+nv)
  PBarkp1 = PBarkp1 + ...
            Wci_p*(xpMat(:,ii) - xBarkp1)*(xpMat(:,ii) - xBarkp1)';
end
ekp1 = xBarkp1(7:9);
RBIBarkp1 = euler2dcm(ekp1)*RBIHatk;                                        % Propagate DCM attitude forward
xBarkp1(7:9) = zeros(3,1);
% Set k = kp1 in preparation for next iteration
RBIBark = RBIBarkp1; xBark = xBarkp1; PBark = PBarkp1;                      % Propagate states forward







