stim_sens = 'stimM2_sensSTN';
% stim_sens = 'stimSTN_sensM2'

for ctype = 1:2
%     load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '.mat'],'powspec_save','intpow','maxpow','burRate','burdur','burAmp','burPPC','burInt','phaseShift','conStren')
    load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '_' stim_sens '_thresholdFitted.mat'],'powspec_save',...
    'intpow','maxpow','burRate','burdur','burAmp','burAmpMid','burPPC',...
    'durStore','ampStore','ppcStore','siStore',...
    'burInt','phaseShift','conStren'); %,'trajStore')
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];

if ctype == 1
    CON = 1;
elseif ctype ==2
    CON = 3;
end
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
        mod = find(ck_1(CON,:)==1);

phaseShift = rad2deg(phaseShift);

    if ctype == 1
        
        cmap = brewermap(25,'RdBu');
        panInd = 1:3;
    else
        cmap = brewermap(25,'RdBu');
        panInd = 4:6;
    end
    cmap(10,:) = [0 0 0]; % fix zero to black
    intpow_S = squeeze(maxpow(2,3,:,:,2:end)); %[M2 STN coh];[lb hb b],[nostim stim],phase,constren
    intpow_nmz = [];
    for K = 1:size(intpow_S,3)
        X = squeeze(intpow_S(:,:,K))';
        %                  X = squeeze(burdur(1,:,:,K))';
        intpow_nmz(:,K) = 100.*(X(:,2)-X(:,1))./X(:,1);
    end
%     intpow_nmz = intpow_nmz(:,2:end);
    
    figure(1)
    NcS_sel = [1:3:19];

    subplot(2,3,panInd(1))
    p = plot(phaseShift,intpow_nmz(:,NcS_sel),'LineWidth',2);
    for i = 1:numel(p)
        p(i).Color = cmap(NcS_sel(i),:);
    end
    grid on
    
    a = gca;
    a.XTick = rad2deg(([0 pi/2 pi 3*pi/2 2*pi]));
    xlim([0 360])
    xlabel('Stimulating Phase')
    ylabel('% Change in Beta Power')
    if CON == 1
    ylim([-100 500]);
    elseif CON == 3
    ylim([-55 300]);
    end        
    title('ARC for 14-30 Hz')
    
    
    subplot(2,3,panInd(2))
    betalist = ck_1(CON,:)*100; %linspace(10,190,19);
%     conStren(2:end)
    p = plotPRCSumStats(betalist,max(intpow_nmz),min(intpow_nmz),range(intpow_nmz),NcS_sel,cmap);
    title('ARC Statistics')
        if CON == 1
    ylim([-100 500]);
    elseif CON == 3
    ylim([-55 300]);
    end     
    
    % Find Phases
    [a Prom] = max(intpow_nmz);
    pp = phaseShift(Prom);
%     pp(pp==2*pi) = 0; % Ensure wt
    [a Sup] = min(intpow_nmz);
    sp = phaseShift(Sup);
%     sp(sp==2*pi) = 0; % Ensure wrapped
    
    subplot(2,3,panInd(3))
    [p sc] = plotPRC_phases(betalist,sp,pp,NcS_sel,cmap)
    title('ARC Phases')
    
    set(gcf,'Position',[470         315        1114         663])
    
    
%     figure(2)
%     plot(R.frqz,squeeze(powspec_save(:,1,1,1,6:end)))
end