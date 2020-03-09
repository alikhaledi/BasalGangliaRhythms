
for ctype = 1:2
    load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '.mat'],'powspec_save','intpow','maxpow','burRate','burdur','burAmp','burPPC','burInt','phaseShift','conStren')
    
    if ctype == 1
        cmap = brewermap(18,'Blues');
        panInd = 1:3;
    else
        cmap = brewermap(18,'Reds');
        panInd = 4:6;
    end
    
    intpow_S = squeeze(maxpow(2,3,:,:,:));
    intpow_nmz = [];
    for K = 1:size(intpow_S,3)
        X = squeeze(intpow_S(:,:,K))';
        %                  X = squeeze(burdur(1,:,:,K))';
        intpow_nmz(:,K) = 100.*(X(:,2)-X(:,1))./X(:,1);
    end
    intpow_nmz = intpow_nmz(1:14,2:end);
    
    NcS_sel = 2:4:18;
    subplot(2,3,panInd(1))
    p = plot(phaseShift,intpow_nmz(:,NcS_sel),'LineWidth',2);
    for i = 1:numel(p)
        p(i).Color = cmap(NcS_sel(i),:);
    end
    grid on
    
    a = gca;
    a.XTick = ([0 pi/2 pi 3*pi/2 2*pi]);
    xlabel('Stimulating Phase')
    ylabel('% Change in Beta Power')
    ylim([-50 125]);
    title('ARC for 14-30 Hz')
    
    
    subplot(2,3,panInd(2))
    p = plotPRCSumStats(conStren(2:end),max(intpow_nmz),min(intpow_nmz),range(intpow_nmz),NcS_sel,cmap);
    title('ARC Statistics')
    
    
    % Find Phases
    [a Prom] = max(intpow_nmz);
    pp = phaseShift(Prom);
    pp(pp==2*pi) = 0; % Ensure wrapped
    
    [a Sup] = min(intpow_nmz);
    sp = phaseShift(Sup);
    sp(sp==2*pi) = 0; % Ensure wrapped
    
    subplot(2,3,panInd(3))
    [p sc] = plotPRC_phases(conStren(2:end),sp,pp,NcS_sel,cmap)
    title('ARC Phases')
    
    set(gcf,'Position',[470         315        1114         663])
end