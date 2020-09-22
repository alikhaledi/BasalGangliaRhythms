function [R,BB] = computeBurstWrapper_V2(R)
% load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')
% R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
R.condname = num2cell(1:20);
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 4 8 16 4 18],:);
hdext = {'','_F1','_bKF','_KrSel'};
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
% R.condcmap(6,:) = [0 0 0];
for CON = [1 3]
    HDM = 3; % This is the scale decided from the 10% 100% and 190% intervals
    BB = [];
    R.condname = {};
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext{HDM} '.mat'],'ck_1')
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim' hdext{HDM} '.mat']); % The low density estimate
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext{HDM} '.mat']); % The low density estimate
    
    if HDM == 1 || HDM == 2 || HDM == 3
        MList = 1:numel(xsim);
    elseif HDM == 4
        MList = R.betaKrange(:,CON)';
    end
    
    xsimMod = {}; bpowr = []; fpow = []; R.condname = {};
    ip = 0;
    for i = MList
        ip = ip+1;
        xsimMod{ip} = 1e7.*xsim{i}{1};
        [bpowr(ip) b] = max(feat{i}(1,4,4,3,:));
        fpow(ip) =R.frqz(b);
        R.condname{ip} = num2str(ck_1(i),2);
    end
    
    
    [R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
    R.BB.thresh_prctile = 75;% o85; tl 80
%     BB.threshold_type = 'localThresh';
        BB.threshold_type = 'baseModelThresh';

    switch BB.threshold_type
        case 'localThresh'
            condsel = MList;
        case 'baseModelThresh'
            condsel = find(ck_1(CON,:)==1); % find base model
    end
    BB = compute_BurstThreshold(R,BB,condsel,0);
    R.BB.minBBlength = 1; %  Minimum burst period- cycles
    BB.plot.durlogflag = 0;
    if HDM == 1 || HDM == 3 || HDM == 4
        memflag = 0;
    else
        memflag = 1; %Crop time series data (for memory)
    end
    BB = defineBetaEvents(R,BB,memflag,0);
    BB = getBurstStatsXConds(BB,0);
    BB.condlist = ck_1;
    save([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) hdext{HDM} BB.threshold_type '.mat'],'BB','-v7.3')
end
end