function [uexs,R,phi] = zeroCrossingPhaseStim_v3(uexs,R,tstep,xstore,dt,uvar,phi)
if (tstep+fix(R.IntP.phaseStim.stimlength/dt))<size(uexs,1) || tstep==0
    if R.IntP.phaseStim.eps == 0
        BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),:);
        
    else
        BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),end-((R.IntP.phaseStim.buff)/dt):end);
    end
    
    if R.IntP.phaseStim.filtflag == 0
        [dum,B,A] = ft_preproc_bandpassfilter(BU, 1/dt, [14 21],4,'but','twopass');
        R.IntP.phaseStim.filtA = A;
        R.IntP.phaseStim.filtB = B;
        R.IntP.phaseStim.filtflag = 1;
    end
    BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,padarray(BU,[0 1/dt]));
    BUB([1:1/dt 1+end-1/dt:end]) = [];
    BEnv = abs(hilbert(BUB));
    
    % BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
    if R.IntP.phaseStim.eps == 0
        R.IntP.phaseStim.eps = prctile(BEnv,R.IntP.phaseStim.epsthresh);
        return
    end
    
    %  ts = 0:dt:3;
    A = (R.IntP.phaseStim.stimAmp*uvar); % setup the amplitude of the stim
    % Check the Envelope for gating
    gate = all(BEnv(end-fix(R.IntP.phaseStim.trackdelay/dt)-fix(R.IntP.phaseStim.minBS/dt):end-fix(R.IntP.phaseStim.trackdelay/dt)) > R.IntP.phaseStim.eps);
    brake = ~all(uexs(tstep-(R.IntP.phaseStim.stimGap/dt):tstep,R.IntP.phaseStim.sensStm(2))==0); % If within
    % This sets up the stim period and ensures breaks
    if gate && ~brake
        uexs(tstep:tstep+(R.IntP.phaseStim.stimPeriod/dt),R.IntP.phaseStim.sensStm(2)) = 1e-32;
    end
    off = uexs(tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) == 0;
    
    % if gate
    % a = 1;
    % end
    if  ~off  %&& gate
        zci = find(diff(sign(BUB))>1,1,'last'); % location of last positive zero crossing
        zstep = (length(BUB)-zci); % location relative to curren time step
        
        % zstep*dt is the current time of the sinusoid relative to t = 0;
        cph = zstep*dt;
        phiPred = 2*pi*18*(cph:dt:cph+R.IntP.phaseStim.stimlength)- pi/2;
        
        stim = sin(phiPred + R.IntP.phaseStim.phaseshift).*A;
        uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = stim;
        phi(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = phiPred;
    end
end
% if brake
%         uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
%     phi(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
% end

% % Demo only
% % These are for demo only
% BA = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(2)),:);
% BAB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,padarray(BA,[0 1/dt]));
% BAB([1:1/dt 1+end-1/dt:end]) = [];
% BAEnv = abs(hilbert(BAB));
% figure(1)
% clf
% plot(BEnv(end-fix(0.3/dt)-fix(R.IntP.phaseStim.minBS/dt):end-fix(R.IntP.phaseStim.trackdelay/dt)))
% hold on
%  plot([0 numel(uexs)],[R.IntP.phaseStim.eps R.IntP.phaseStim.eps],'k--')
%  xlim([0 600])
%  drawnow
% a = 1;
% figure(2)
% clf
% yyaxis left; plot(BAEnv); hold on;plot(BA); plot(BAB); yyaxis right; plot(uexs(:,R.IntP.phaseStim.sensStm(2)))
% yyaxis left; plot([0 numel(uexs)],[R.IntP.phaseStim.eps R.IntP.phaseStim.eps],'k--')
% xlim([tstep-2e3 tstep+2e3])
% ylim([-1.5 1.5]*1e-6);
% drawnow