function [uexs,eps,R] = zeroCrossingPhaseStim_v2(uexs,R,tstep,xstore,dt,eps,uvar)
if eps == 0
    BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),:);
    
else
    BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),end-25:end);
end
BU = BU-mean(BU);
% BUB = bandpass(BU,[15 21],1/dt);
% BUB2  = filtfilt(R.IntP.phaseStim.bpfilt,BU);
if R.IntP.phaseStim.filtflag == 0
    %     [dum,B,A] = ft_preproc_bandpassfilter(BU, 1/dt, [14
    %     21],4,'but','twopass');
    %     R.IntP.phaseStim.filtA = A;
    %     R.IntP.phaseStim.filtB = B;
    
    SR = 1/dt;
    BW = 2;
    wo = 18/(SR/2);
    bw = BW/(SR/2);
    qf =3;  % Q factor ( -qf);
    
    % [bnc,anc] = iirnotch(wo,bw,qf);
    [bpk,apk] = iirpeak(wo,bw,qf);
    
    R.IntP.phaseStim.filtA = apk;
    R.IntP.phaseStim.filtB = bpk;
    R.IntP.phaseStim.filtflag = 1;
end
BUB = filter(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,BU);

BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
if eps == 0
    eps = prctile(BEnv,10);
    return
end
% These are for demo only
%     BA = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(2)),:);
% BAB = filter(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,BA);
% BAEnv = smooth(abs(BAB),200);


%  ts = 0:dt:3; 
if all(BEnv(end-50-fix(R.IntP.phaseStim.minBS/dt):end-50) > eps)
%     if uexs(tstep,R.IntP.phaseStim.sensStm(2)) ==0
        
        zci = find(diff(sign(BUB))>1,1,'last'); % location of last positive zero crossing
        zstep = (length(BUB)-zci); % location relative to curren time step
        
        % zstep*dt is the current time of the sinusoid relative to t = 0;
        cph = zstep*dt;
        phiPred = 2*pi*18*(cph:dt:cph+R.IntP.phaseStim.stimlength);
        
        stim = sin(phiPred + R.IntP.phaseStim.phaseshift).*(R.IntP.phaseStim.stimAmp*uvar);
        uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = stim;
%     else
%         uexs = uexs;
%     end
else
    uexs(tstep+1:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
    
%     a = 1;
end
% Demo only
% clf
% plot(BEnv(1:end-50)); drawnow; %BEnv(1:end-fix(R.IntP.phaseStim.trackdelay/dt))); drawnow
% drawnow
% a = 1;
% clf
% yyaxis left; plot(BAEnv); hold on;plot(BA); plot(BAB); yyaxis right; plot(uexs(:,R.IntP.phaseStim.sensStm(2)))
% xlim([tstep-2e3 tstep+2e3])
% drawnow