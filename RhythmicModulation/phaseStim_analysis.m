clear all;
load('tmp_sesh')
load([R.rootn 'data\rat_InDirect_ModelComp\phaseStimSave\stim_tmp'],'uexs')
pU = uexs(4,round(R.obs.brn*(1/R.IntP.dt))+1:end);
betaBurstInds = SplitVec(find(abs(pU)>0),'consecutive'); % Split up data based upon the target threshold
segL = cellfun('length',betaBurstInds);
% Now crop using minimum
fsamp = (1/R.IntP.dt);
burstSelection = find(segL>((1.5/18)*fsamp));


%% Get the data features (bandpass
X = xsim_ip{2}{1}([1 4],:);
XF = ft_preproc_bandpassfilter(X,fsamp,[14 21],[],'fir');
XEnv = abs(hilbert(XF));
XPhi = angle(hilbert(XF));

winsize = 0.3*fsamp;

for seg = 1:numel(burstSelection)
    zeroP = betaBurstInds{burstSelection(seg)}(1);
    befWin = zeroP-winsize:zeroP;
    aftWin = zeroP:zeroP+winsize;

    befEnv(:,seg) = XEnv(2,befWin);
    aftEnv(:,seg) = XEnv(2,aftWin);
    
end
