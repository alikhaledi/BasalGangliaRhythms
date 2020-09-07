function [uexs,eps,R] = phaseRegressionStim(uexs,R,tstep,xstore,dt,eps,uvar)


%                 prctile(xstore(R.obs.outstates(4),R.obs.brn/dt:end);
% BUL = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),:); %end-((R.IntP.phaseStim.buff)./dt):end);
% BULB = BUL-mean(BUL);
% BULB  = filtfilt(R.IntP.phaseStim.filtb,R.IntP.phaseStim.filta,BULB);
% %                 %                 BUB = bandpass(BU,[18 24],1/dt);
% %                 HB = hilbert(BUB);
% % BULEnv = smooth(abs(BULB),100);
% BULBA = abs(hilbert(BULB));

%% Starts here
BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm),end-((R.IntP.phaseStim.buff)/dt):end);
BU = BU-mean(BU,2);
% Initialize filters
if R.IntP.phaseStim.filtflag == 0
    [BUB,B,A] = ft_preproc_bandpassfilter(BU, 1/dt, [15 21]);
    R.IntP.phaseStim.filta = A;
    R.IntP.phaseStim.filtb = B;
    R.IntP.phaseStim.filtflag = 1;
end

% Start with buffered data
BUB = []; BEnv = [];
for i = 1:2
    BUB(i,:)  = filtfilt(R.IntP.phaseStim.filtb,R.IntP.phaseStim.filta,BU(i,:));
end
BPhi = angle(hilbert(BUB(2,:))); % phase of stim pop
BEnv = smooth(abs((BUB(1,:))),800); %abs(HB); %smooth(abs(HB),100); % env of sens

if eps == 0
    eps = prctile(BEnv,75);
end

if any(BEnv(end-fix(R.IntP.phaseStim.minBS/dt):end) > eps)
    %     if uexs(tstep,2) ==0
    % Perform local regression on phase
    Y =  unwrap(BPhi(end-fix(R.IntP.phaseStim.regleng/dt):end))'; % Phase of stimulated signal
    X = (tstep-fix(R.IntP.phaseStim.regleng/dt):1:tstep)';
    X = [ones(length(X),1) X];
    b = X\Y;
    %                                         yCalc = X*b;
    X2 = (tstep+1:1:tstep+fix(R.IntP.phaseStim.stimlength/dt))';
    X2 = [ones(length(X2),1) X2];
    phiPred = X2*b;
    
    stim = wrapToPi(phiPred + R.IntP.phaseStim.phaseshift); %sin(phiPred + R.IntP.phaseStim.phaseshift).*(R.IntP.phaseStim.stimAmp*uvar);
    uexs(tstep+1:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = stim;
    %     else
    %         uexs = uexs;
    %     end
else
    uexs(tstep+1:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
    
    a = 1;
end


%                  yyaxis left; plot(BULEnv); yyaxis right; plot(uexs(:,4))
%                  drawnow