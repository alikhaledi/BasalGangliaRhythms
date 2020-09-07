clear; close all

R.IntP.phaseStim.sensStm = [1 1];
R.IntP.dt = 5e-4;
R.IntP.buffer = 50;
R.obs.brn = 3;
R.obs.outstates = 1;


R.IntP.phaseStim.switch = 1;
%                 R.IntP.phaseStim.stimfreq = 18;
R.IntP.phaseStim.buff = 3; % This is the buffer used to compute the current phase
R.IntP.phaseStim.minBS =  ((1/18)*(1./R.IntP.dt))/1000; % Minimum burst length
R.IntP.phaseStim.trackdelay = 0.25; % this is the delay to take (as the end of the hilbert in unstable
R.IntP.phaseStim.stimlength = 0.3; % 300ms stim delivery
R.IntP.phaseStim.stimAmp = 1/4; % times the variance of the normal input
R.IntP.phaseStim.regleng = 3/18; % 500ms regression only
%                 R.IntP.phaseStim.thresh = BB.epsAmp;
R.IntP.phaseStim.bpfilt = designfilt('bandpassiir', 'FilterOrder', 20,...
    'HalfPowerFrequency1', 15, 'HalfPowerFrequency2', 21,...
    'SampleRate', 1/R.IntP.dt);
R.IntP.phaseStim.phaseshift = 0;
R.IntP.phaseStim.filtflag = 0;

%
R.IntP.phaseStim.minBS = 0;
R.IntP.phaseStim.epsthresh = 0;

%% Start some tests

for test =1:2
    
    if test == 2
        load exampledata
    elseif test == 1
        tData=80;  % in Secs
        SR=1/R.IntP.dt;
        nData=SR*tData + 1;
        tmaxis=(0:nData-1)*R.IntP.dt;
        
        fr1=18;
        am1=1;
        am2=2;
        
        dt1=am1*sin(2*pi*fr1*tmaxis+0);
        dt2=am2*(rand(1,nData)-0.5);
        dta=dt1+dt2;
        tmpdata = dta;
    end
    for res = 1:2
        if res == 1
            rl = 5;
        elseif res == 2
            rl = 25;
        end
        epsStim = zeros(numel(tmpdata),1); dt = R.IntP.dt; uexs = zeros(numel(tmpdata),1); phi = zeros(numel(tmpdata),1);
        for tstep = R.IntP.buffer:numel(tmpdata)
            X = tmpdata(1:tstep);
            
            
            if tstep >((R.obs.brn)/dt) && (rem(tstep,rl) == 0) %&& ~any(uexs(tstep,:))
                if R.IntP.phaseStim.switch
                    [uexs,epsStim,R,phi] = zeroCrossingPhaseStim(uexs,R,tstep,X,dt,epsStim,std(tmpdata),phi);
                    %                 [uexs,epsStim,R] = phaseRegressionStim(uexs,R,tstep,xstore,dt,epsStim,std(us(:,R.IntP.phaseStim.sensStm(2))));
                end
                sprintf('%0.3f',tstep/numel(tmpdata))
            end
        end
        
        % Get true phase
        BU = tmpdata;
        BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,padarray(BU,[0 1/dt]));
        BUB([1:1/dt 1+end-1/dt:end]) = [];
        
        phiT = angle(hilbert(BUB));
        
        pred = wrapToPi(phi(1:numel(phiT)))';
        act = phiT;
        
        PLV(test,res) = abs(sum(exp(i*(pred-act))))./numel(act);
        tvec = linspace(0,numel(pred)*dt,numel(pred));
        subplot(2,2,sub2ind([2 2],res,test)); plot(tvec,pred); hold on; plot(tvec,act)
        xlim([75 76])
    end
end

figure
b = bar(PLV);
a = gca;
a.XTickLabel = {'Noisy Sine','Sim. Network Activity'}
ylabel('Phase Locking Value between predicted-actual')
legend({'5ms refresh','25ms refresh'})
