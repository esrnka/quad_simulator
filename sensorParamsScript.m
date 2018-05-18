% sensorParamsScript.m
%
% Loads quadrotor sensor parameters into the structure sensorParams

%-------------- GNSS ---------------
% rA(:,1) holds the 3x1 position of the primary GNSS antenna in the body frame,
% in meters.  rA(:,2) holds the position of the secondary antenna.
sensorParams.rA = [0.1013 -0.0887; -0.0004 -0.0004; 0.0472 0.0472];
% The 3x3 error covariance matrix for the unconstrained precise GNSS (RTK)
% solution for either the primary or secondary antenna, in ENU coordinates, in
% meters^2.
sensorParams.RpL = diag([0.006^2, 0.006^2, 0.012^2]);
% Standard deviation of the error in the measured baseline-constrained vector
% from the primary to the secondary antenna, in meters.
sensorParams.sigmaC = 0.012;
% 3x1 position of the origin of the I frame, which is a local ENU frame,
% expressed in ECEF coordinates in meters
sensorParams.r0G = [-742015.136; -5462219.538; 3198014.35];
% 3x1 position of reference antenna expressed in the ECEF frame, in meters
sensorParams.rRefG = [-741990.536; -5462227.638; 3198019.45];

%----------- HD Camera -------------
% 3x1 position of the HD camera center in the body frame, in meters.
sensorParams.rc = [0.1159; -0.0004; -0.0435];
% 3x3 direction cosine matrix relating the HD camera and body reference
% frames.  For some vector v, vC = RCB*vB.
sensorParams.RCB = euler2dcm([pi/2, 120*pi/180, 0]);
% 2x2 error covariance matrix for the Gaussian image noise, in pixels^2
sensorParams.Rc = diag([20^2,20^2]);
% Distance of the image plane along camera z axis, in meters
sensorParams.f = 0.004;
% Camera intrinsic matrix
sensorParams.K = diag([sensorParams.f,sensorParams.f,1]);
% Pixel size, in meters
sensorParams.pixelSize = 2e-6;
% Image plane size (i.e., width and height), in meters, along the camera x and
% y dimensions
sensorParams.imagePlaneSize = [0.004,0.004];

%------------- IMU -----------------
% IMU sampling interval, in seconds
sensorParams.IMUdelt = 0.005;
% 3x1 position of IMU accelerometer proof mass in the B frame, in meters.
sensorParams.lB = [-0.0884; 0.0134; -0.0399];
% Accelerometer white noise from accel. spec. sheet, in (milli-g)^2/Hz
sensorParams.Sa = (0.1)^2; 
% Error covariance matrix for accelerometer noise, in (m/s^2)^2
sensorParams.Qa = (9.8/1000)^2*(sensorParams.Sa/sensorParams.IMUdelt)*eye(3);
% Accelerometer bias time constant, in seconds
sensorParams.taua = 100;  
% Accelerometer bias from accel. spec. sheet, in milli-g
sensorParams.sigmaa = 10; 
% Error covariance matrix for accelerometer bias noise, in (m/s^2)^2
sensorParams.alphaa = exp(-sensorParams.IMUdelt/sensorParams.taua); 
sensorParams.Qa2 = (9.8/1000)^2*(sensorParams.sigmaa)^2*(1 - sensorParams.alphaa^2)*eye(3);
% Gyro white noise from spec. sheet, in (deg/s)^2/Hz
sensorParams.Sg = (1e-2)^2; 
% Error covariance matrix for gyro noise, in (rad/s)^2
sensorParams.Qg = (pi/180)^2 * (sensorParams.Sg/sensorParams.IMUdelt)*eye(3);
% Gyro bias time constant, in seconds
sensorParams.taug = 100; 
% Gyro bias from spec. sheet, in (deg/h)
sensorParams.sigmag = 10; 
% Error covariance matrix for gyro bias noise, in (rad/s)^2
sensorParams.alphag = exp(-sensorParams.IMUdelt/sensorParams.taug);
sensorParams.Qg2 = (pi/(3600*180))^2*(sensorParams.sigmag)^2*(1-sensorParams.alphag^2)*eye(3);
