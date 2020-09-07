function R = cohMechanism(R)
close all
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];

% load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat_F1.mat'],'feat_HD','feat_STR_GPe')
R.CONnames = {'M2 -> STN','STR -| GPe','GPe -| STN','STN -> GPe'};
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
cmap1 = brewermap(40,'Reds');
% cmap1(22,:) = [0 0 0];
cmap2 = brewermap(40,'Blues');

ip = 0;
for CON = [1 3]
        ip = ip + 1;
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim_F1.mat'])

    
    
end