function [R] = BAA_plotStimTypeOutcomes_Eps(Rorg)
state = 1; % baseline state only
CON = 1; % dont need to discriminate
close all
SStype = {'stimM2_sensSTN','stimSTN_sensM2','stimSTN_sensGPe','stimM2_sensSTN','stimSTN_sensM2'};
phaseN = [12];
EpsVec = [0:15:100];

cmap{1} = [linspace(1,0,7); zeros(1,7); zeros(1,7)]';
cmap{2} = [zeros(1,7); zeros(1,7); linspace(1,0,7)]';

for EpsN = 1:7
    SScomb = 1
    
    %% Loop through Connections
    rootan = [Rorg.rootn 'data\phaseLockedStimEpsComp'];
    load([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_feat' num2str(SScomb) '.mat'],'feat_sim_save')
    load([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_xsim' num2str(SScomb) '.mat'],'xsim_ip')
    load([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_Rout' num2str(SScomb) '.mat'],'Rout')
    load([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_pU_save' num2str(SScomb) '.mat'],'pU_save')
    
    baseFeat = squeeze(feat_sim_save{1,EpsN}{1}(1,4,4,1,:));
    bp_base(1) = bpSpec(Rout.frqz,baseFeat,[14 21]);
    bp_base(2) = bpSpec(Rout.frqz,baseFeat,[21 30]);
    state = 1; %:size(ck_1,2)
    deltaBP = []; phaseFeat = [];
    for p = 1:phaseN(SScomb)
        phaseFeat(:,p) = squeeze(feat_sim_save{2,EpsN}{p}(1,4,4,1,:));
        bp = bpSpec(Rout.frqz,phaseFeat(:,p),[14 21]);
        deltaBP(1,p) = 100*(bp-bp_base(1))./bp_base(1);
                bp = bpSpec(Rout.frqz,phaseFeat(:,p),[21 30]);

                deltaBP(2,p) = 100*(bp-bp_base(2))./bp_base(2);

    end
    [~,supp] = min(deltaBP(1,:));
    [~,amp] = max(deltaBP(1,:));
    subplot(1,3,1)
    p1(1) = plot(Rout.frqz,baseFeat,'g'); hold on
    p1(EpsN+1) = plot(Rout.frqz,phaseFeat(:,supp),'Color',cmap{1}(EpsN,:))
       xlabel('Hz'); ylabel('Power'); xlim([2 48]); box off; axis square
       
    subplot(1,3,2)
    p2(1) =  plot(Rout.frqz,baseFeat,'g'); hold on
    p2(EpsN+1) = plot(Rout.frqz,phaseFeat(:,amp),'Color',cmap{2}(EpsN,:))
    deltaBPStore(:,:,EpsN) = deltaBP(:,[supp amp])
    
    xlabel('Hz'); ylabel('Power'); xlim([2 48]); box off; axis square
    
    
end

subplot(1,3,1)
legend(p1([1 3 5 7]) ,{'Base','15th','45th','75th %tile'})
title('Suppressing')
subplot(1,3,2)
legend(p2([1 3 5 7]) ,{'Base','15th','45th','75th %tile'})
title('Amplifying')

subplot(1,3,3)
plot(EpsVec,squeeze(deltaBPStore(1,1,:)),'r')
hold on
plot(EpsVec,squeeze(deltaBPStore(1,2,:)),'b')
plot(EpsVec,squeeze(deltaBPStore(2,1,:)),'r--')
hold on
plot(EpsVec,squeeze(deltaBPStore(2,2,:)),'b--')
title('Summary')
    xlabel('Threshold (percentile)'); ylabel('Beta BP Change %'); box off; axis square
    legend({'LB Pow (amp)','LB Pow (sup)','HB Pow (amp)','HB Pow (sup)'})
set(gcf,'Position',[ 290         235        1490         743])
a = 1;

% ! shutdown /s

function bp = bpSpec(Hz,Pxy,flim)
bandInds = find(Hz>=flim(1) & Hz<= flim(2));
bp = sum(Pxy(bandInds))*(numel(bandInds)*(Hz(2)-Hz(1)));