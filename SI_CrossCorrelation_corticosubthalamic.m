function R = SI_CrossCorrelation_corticosubthalamic(R)
close all
rootan = [R.rootn 'data\ConnectionSweep'];

R.CONnames = {'M2 -> STN','GPe -| STN'};
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
cmap1 = brewermap(40,'Greys');
cmap2 = brewermap(40,'Blues');
scmap = brewermap(4,'Set1');
statecmap{1} = scmap(1:2,:);
statecmap{2} = scmap(3:4,:);
hdext = '_REV'; % version tag

for CON = 2
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext '.mat'])
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext '.mat'])
    
    
    
    figure(1)

end

