% Overlap subplots
close all
    load([rootan '\BB_' R.out.tag '_stateConnectivityMatrix.mat'],'connectMat','specMat','statBurstOverl','sampBurstOverl')

figure(101)
ha =   tight_subplot(3,3,0.05);
delete(ha([2 4 6 8]))
R.statename = {'Fitted','HD-Down','HD-Up';'' '' ''; 'Fitted','PS-Down','PS-Up'};

scmap = brewermap(4,'Set1');
for CON = [1 3]
        if CON == 1
        cmap = brewermap(128,'RdBu');
    elseif CON == 3
        cmap = brewermap(128,'PRGn');
        end
    
        figure(101)
        plotOverlapMats(R,statBurstOverl,CON,cmap,ha)
end

set(gcf,'Position',[652   135   860   646])

figure(102)
plotOverlapBars(sampBurstOverl)