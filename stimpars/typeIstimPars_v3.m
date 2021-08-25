function R = typeIstimPars_v3(R)
R.IntP.phaseStim.filtflag = 0;
R.IntP.phaseStim.buff = 1; % This is the buffer used to compute the current phase
R.IntP.phaseStim.minBS = 5/1000; % WAS 0 Minimum burst length
R.IntP.phaseStim.trackdelay = 25/1000; % WAS 25; % this is the delay to take (as the end of the hilbert in unstable
R.IntP.phaseStim.upperiod = 25/1000; % WAS 10ms update period
R.IntP.phaseStim.stimlength = R.IntP.phaseStim.upperiod; %5/1000; % 300ms stim delivery
R.IntP.phaseStim.stimAmp = 1/3; %WAS 1/3; % times the variance of the normal input; % Might need higher for STN stim
R.IntP.phaseStim.stimPeriod = 0.5;  %WAS 0.5 % stimulation period
R.IntP.phaseStim.stimGap = 1; % break between stimulation

R.IntP.phaseStim.filtflag = 0;
R.IntP.phaseStim.epsthresh = 75; 
R.IntP.phaseStim.eps = 0;