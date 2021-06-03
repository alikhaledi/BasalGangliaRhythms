function [burstSelInds,XF,XEnv,XPhi,epsAmp] = burstOverlap(X,fsamp,frqz,minper)
% Filter Data
[dum,bpk,apk] = ft_preproc_bandpassfilter(X, fsamp,frqz,4,'but','twopass');
XF =filtfilt(bpk,apk,X')';


% Analytic Signal
XEnv = []; XPhi = [];
for ch = 1:6
    XEnv(ch,:) = abs(hilbert(XF(ch,:)));
    XPhi(ch,:) = angle(hilbert(XF(ch,:)));
    % Set Threshold
    epsAmp(ch) = prctile(XEnv(ch,:),75);
end

% Define Bursts
ThreshX = double(XEnv(4,:) > epsAmp(4));
minS = (minper/frqz(1))*fsamp; % set minimum to 3 periods of slowest
betaBurstInds = SplitVec(find(ThreshX),'consecutive'); % Split up data based upon the target threshold
segL = cellfun('length',betaBurstInds); % Find burst lengths
burstSelInds = segL>minS; % Select bursts with above min length
burstSelInds = betaBurstInds(burstSelInds);
% 
