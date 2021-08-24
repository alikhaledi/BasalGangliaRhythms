function R = plotSweepSpectraBasic(R)
close all
rootan = [R.rootn 'data\' R.out.tag '\ConnectionSweep'];

R.CONnames = {'M2 -> STN','GPe -| STN'};
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
cmap1 = brewermap(40,'Greys');
cmap2 = brewermap(40,'Blues');
scmap = brewermap(4,'Set1');
statecmap{1} = scmap(1:2,:);
statecmap{2} = scmap(3:4,:);
hdext = '_REV'; % version tag

for CON = 1:2
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext '.mat'])
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext '.mat'])
    figure(1)
    if CON == 1
        subplot(2,3,1)
    elseif CON == 2
        subplot(2,3,4)
    end
    
    a = plotSweepSpectra(R.frqz,feat,feat{1},cmap1,{R.condname{[2 1 3]}}, [1 15 30],2:4:31,[1,1,1],statecmap{CON})
    title('M2 power')
        ylim([1e-15 5e-14])
        set(gca, 'YScale', 'log'); %, 'XScale', 'log')
    if CON == 1
        subplot(2,3,2)
    elseif CON == 3
        subplot(2,3,5)
    end
    plotSweepSpectra(R.frqz,feat,feat{1},cmap1,{R.condname{[2 1 3]}}, [1 15 30],2:4:31,[4,4,1],statecmap{CON})
    title('STN power')
        ylim([1e-16 1e-12])
        set(gca, 'YScale', 'log'); %, 'XScale', 'log')
    
    if CON == 1
        subplot(2,3,3)
    elseif CON == 3
        subplot(2,3,6)
    end
    plotSweepSpectra(R.frqz,feat,feat{1},cmap1,{R.condname{[2 1 3]}}, [1 15 30],2:4:31,[4,1,4],statecmap{CON})
    title('M2/STN coherence')
            ylim([0 1])
    
end

for CON = 1:2
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext '.mat'])
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext '.mat'])
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim' hdext '.mat'])
    bpow = []; fpow = [];
    [bpowr_br,fpow_br,bpowr,fpow,bcohr,fcoh,fpowCtx,bpowrCtx] = computeBetaSpectralStats(R.frqz,feat);
    ck_1 = ck_1(CON,:); % The scale for this connection modification
    ck_1 = ck_1(2:end);

    % Scale bpow to 0 (for plotting)
    [a zind] = min(abs(ck_1-1)); % base model
    bpowr = 100*(bpowr(2:end)-bpowr(1))/bpowr(1);
    bpowrCtx = 100*(bpowrCtx(2:end)-bpowrCtx(1))/bpowrCtx(1);
    bcohr = 100*(bcohr(2:end)-bcohr(1))/bcohr(1);
    fpow = fpow(2:end);
    fpowCtx = fpowCtx(2:end);
    fcoh = fcoh(2:end);
    
    powInds = find(bpowr>500);
    fpowCtx(powInds) = nan(1,numel(powInds));
    bpowrCtx(powInds) = nan(1,numel(powInds));
    fcoh(powInds) = nan(1,numel(powInds));
    bcohr(powInds) = nan(1,numel(powInds));
    
    Krange{CON} = ck_1(bpowr<500);
    
    bsel = 2:4:31;
    
    %% Data Selection
    indsel = [bsel(1) bsel(end)];
    dataSelect{CON} = xsim([1 bsel(1) bsel(end)]);
    dataProperties(:,:,CON) = [bpowrCtx(indsel); fpowCtx(indsel); bpowr(indsel); fpow(indsel);  bcohr(indsel); fcoh(indsel);] 
    
end
save([rootan '\BB_' R.out.tag '_DiscreteData.mat'],'dataSelect','dataProperties')

% legend(Sr,ck_1(1,:)*100)
R.Krange = Krange;
save([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_KRange.mat'],'Krange')
% R.betaKrange(3,3) = 19;
set(gcf,'Position',[ 518         250        1211         633])