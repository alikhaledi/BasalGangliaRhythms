function BB = compute_BurstThreshold(R,BB,condsel,plotop)
surflag = 0;
% Find the amplitude threshold from the prctile of concatanated data
switch R.BB.threshold_type
    case 'concatSimType' % Concatanate per simulation type
        if surflag == 0
            concatA = [BB.AEnv{condsel}];
            BB.epsAmpfull = prctile(concatA(:,:),R.BB.thresh_prctile ,2); % This is for all channels
            BB.epsAmp = prctile(concatA(R.BB.pairInd(2),:),R.BB.thresh_prctile ,2); % Or only the target channel
            concatA = [BB.PLV{condsel}];
            BB.epsPLV = prctile(concatA(1,:),R.BB.thresh_prctile ,2);
        else
            BORG = load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_org'],'BB');
            BB.epsAmp = BORG.BB.epsAmp;
            BB.epsPLV = BORG.BB.epsPLV ;
            clear BORG
        end
    case 'baseModelThresh'
        A = BB.AEnv{condsel};
        BB.epsAmpfull(:,1) = prctile(A,R.BB.thresh_prctile ,2); % This is for all channels
        BB.epsAmp(1) = prctile(A(R.BB.pairInd(2),:),R.BB.thresh_prctile ,2); % This is for all channels
        A = BB.PLV{condsel};
        BB.epsPLV(1) = prctile(A(1,:),R.BB.thresh_prctile ,2);
        
    case 'localThresh'
        for c = condsel
            A = BB.AEnv{c};
            BB.epsAmpfull(:,c) = prctile(A,R.BB.thresh_prctile ,2); % This is for all channels
            BB.epsAmp(c) = prctile(A(R.BB.pairInd(2),:),R.BB.thresh_prctile ,2); % This is for all channels
            A = BB.PLV{c};
            BB.epsPLV(c) = prctile(A(1,:),R.BB.thresh_prctile ,2);
        end
end

% Plot Bursting
if plotop == 1
    figure(1)
    plotExampleBurstPLV(R,BB)
    set(gcf,'Position',[680  358  1048  622])
end

