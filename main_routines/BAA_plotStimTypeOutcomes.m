function [R] = BAA_plotStimTypeOutcomes(Rorg)
state = 1; % baseline state only
     CON = 1; % dont need to discriminate
close all
SStype = {'stimM2_sensSTN','stimSTN_sensM2','stimSTN_sensGPe','stimM2_sensSTN','stimSTN_sensM2'};
phaseN = [12 12 1 1 1 1 1 1 12 12];
% phaseN = [1 1 1 1 1 1 1 1 1];


for SScomb = 1:10
    subplot(5,2,SScomb)
    %% Loop through Connections
        rootan = [Rorg.rootn 'data\phaseLockedStim'];
        load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],'feat_sim_save')
        load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip')
        load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_Rout' num2str(SScomb) '.mat'],'Rout')
        load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_pU_save' num2str(SScomb) '.mat'],'pU_save')

        baseFeat = squeeze(feat_sim_save{1,state}{1}(1,4,4,1,:));
        plot(Rout.frqz,baseFeat,'k'); hold on
        bp_base = bpSpec(Rout.frqz,baseFeat,[14 21]);
        
        state = 1; %:size(ck_1,2)
        deltaBP = []; phaseFeat = [];
        for p = 1:phaseN(SScomb)
            phaseFeat(:,p) = squeeze(feat_sim_save{2,state}{p}(1,4,4,1,:));
            bp = bpSpec(Rout.frqz,phaseFeat(:,p),[14 21]);
            deltaBP(p) = 100*(bp-bp_base)./bp_base;
        end
        if numel(deltaBP)>1
            [~,supp] = min(deltaBP);
                        [~,amp] = max(deltaBP);
            plot(Rout.frqz,phaseFeat(:,supp),'b')
            plot(Rout.frqz,phaseFeat(:,amp),'r')
            legend({'0' num2str(deltaBP(supp)) num2str(deltaBP(amp))})
        else
            plot(Rout.frqz,phaseFeat,'g')
            legend({'0' num2str(deltaBP(1)) })
        end
        xlabel('Hz'); ylabel('Power'); xlim([2 48])
        % put in titles
        if SScomb == 1
            title('M2 Stim / STN Sense')
            ylabel('Phase locked 14-21 Hz')
        elseif SScomb == 2
            title('M2 Stim / STN Sense')
        elseif SScomb == 3
            ylabel('Raw 14-21 Hz')
        elseif SScomb == 5
            ylabel('Pulse 130 Hz')
        elseif SScomb == 7
            ylabel('Sine 130 Hz')
        elseif SScomb == 9
            ylabel('Playback 14-21 Hz')   
        end

% %                set(gca, 'YScale', 'log');
end

a = 1;

% ! shutdown /s

function bp = bpSpec(Hz,Pxy,flim)
bandInds = find(Hz>=flim(1) & Hz<= flim(2));
bp = sum(Pxy(bandInds))*(numel(bandInds)*(Hz(2)-Hz(1)));