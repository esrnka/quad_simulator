clear all; clc;
% Create the empty data matrices
data = []; inputs = []; targets = [];
for i = 1:20
    % Load the appropriate data file
    datafile = ['run'  num2str(i)  '.mat'];
    load(datafile);
    
    % resample the oversampled state data
    rMat = Q.state.rMat(1:2:end-1,:);
    vMat = Q.state.vMat(1:2:end-1,:);
    omegaBMat = Q.state.omegaBMat(1:2:end-1,:);
    eMat = Q.state.eMat(1:2:end-1,:);
    
    % Number of input/target pairs
    N = Q.tVec(end);
    
    % Build the interleaved raw data matrix
    for ii = 1:length(Q.meas.tk)
        data = [data; Q.meas.rPtilde(ii,:)'; Q.meas.rStilde(ii,:)';...
            Q.meas.rCtilde(ii,:)'; Q.meas.fB(ii,:)'; ...
            Q.meas.omegaBtilde(ii,:)'; rMat(ii,:)'; vMat(ii,:)'; ...
            omegaBMat(ii,:)'; eMat(ii,:)'];
    end
    
    % Size of a single input vector
    s = length(data) / N;
    
    % Split the raw data matrix into 1 second samples
    for jj = 0:N-1
        inputs = [inputs, data(s*jj+1:s*(jj+1))];
        targets = [targets, P.sensorParams.lB];
    end
    
    % Cleanup for the next loop
    data = [];
end
save('nndata.mat','inputs','targets');