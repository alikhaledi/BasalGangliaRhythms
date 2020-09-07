function [R,BB] = checkSimulations(R)
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
          bip=       find(ck_1(CON,:)==1)
X = squeeze(feat{bip}(1,4,4,3,:));
plot(R.frqz,X)
hold on
end