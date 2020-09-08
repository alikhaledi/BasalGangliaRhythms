function [uexs,R,phi] = zeroCrossingPhaseStim(uexs,R,tstep,xstore,dt,uvar,phi)
if R.IntP.phaseStim.eps == 0
    BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),:);
    
else
    BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),end-((R.IntP.phaseStim.buff)/dt):end);
end
BU = BU-mean(BU);

if R.IntP.phaseStim.filtflag == 0
    [dum,B,A] = ft_preproc_bandpassfilter(BU, 1/dt, [14 21],4,'but','twopass');
    R.IntP.phaseStim.filtA = A;
    R.IntP.phaseStim.filtB = B;
    R.IntP.phaseStim.filtflag = 1;
end
BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,padarray(BU,[0 1/dt]));
BUB([1:1/dt 1+end-1/dt:end]) = [];

BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
if R.IntP.phaseStim.eps == 0
    R.IntP.phaseStim.eps = prctile(BEnv,R.IntP.phaseStim.epsthresh);
    return
end



%  ts = 0:dt:3;
if all(BEnv(end-50-fix(R.IntP.phaseStim.minBS/dt):end-50) > R.IntP.phaseStim.eps)
    %     if uexs(tstep,R.IntP.phaseStim.sensStm(2)) ==0
    
    zci = find(diff(sign(BUB))>1,1,'last'); % location of last positive zero crossing
    zstep = (length(BUB)-zci); % location relative to curren time step
    
    % zstep*dt is the current time of the sinusoid relative to t = 0;
    cph = zstep*dt;
    phiPred = 2*pi*18*(cph:dt:cph+R.IntP.phaseStim.stimlength)- pi/2;
    
    stim = sin(phiPred + R.IntP.phaseStim.phaseshift).*(R.IntP.phaseStim.stimAmp*uvar);
    uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = stim;
    phi(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = phiPred;
    %     else
    %         uexs = uexs;
    %     end
else
    uexs(tstep+1:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
    
    %     a = 1;
end
% Demo only
% These are for demo only
%     BA = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(2)),:);
% BAB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,BA);
% BAEnv = smooth(abs(BAB),200);

% clf
% plot(BEnv(1:end-50)); drawnow; %BEnv(1:end-fix(R.IntP.phaseStim.trackdelay/dt))); drawnow
% drawnow
% a = 1;
% clf
% yyaxis left; plot(BAEnv); hold on;plot(BA); plot(BAB); yyaxis right; plot(uexs(:,R.IntP.phaseStim.sensStm(2)))
% xlim([tstep-2e3 tstep+2e3])
% drawnow