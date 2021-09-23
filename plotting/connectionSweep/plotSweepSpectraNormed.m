function a=  plotSweepSpectra(Hz,feat,featemp,cmap,legn,legsel,condsel,chsel,stcmap,featNoise)
if nargin<8
    chsel = 4;
end
if nargin<10
    noiseFlag = 0;
else
        noiseFlag = 1;
end
% Mark out frequency borders
plot([14 14],[1e-32 1],'k--')
hold on
plot([21 21],[1e-32 1],'k--')
plot([30 30],[1e-32 1],'k--')


j = 0;
for i = condsel
    j = j+1;
    if max(squeeze(feat{i}(1,4,4,1,:)))<1e-4; % If STN is over reasonable level
        X = squeeze(feat{i}(1,chsel(1),chsel(2),chsel(3),:));
        if noiseFlag == 1
            XN = squeeze(featNoise(1,chsel(1),chsel(2),chsel(3),:));
            X = 10.*log10(X./XN); % normalize
        end
        
        a(i) = plot(Hz,X,'color',cmap(i,:),'LineWidth',2);
        %     a(i) = plot(Hz,squeeze(feat{i}(1,4,1,2,:)),'color',cmap(i,:),'LineWidth',2);
        hold on
        pind(j) = i;
    end
end
a(i+1) = plot(Hz,squeeze(featemp(1,chsel(1),chsel(2),chsel(3),:)),'k','LineWidth',2);
a(i+2) = plot(Hz,squeeze(feat{2}(1,chsel(1),chsel(2),chsel(3),:)),'Color',stcmap(1,:),'LineWidth',2);
a(i+3) = plot(Hz,squeeze(feat{31}(1,chsel(1),chsel(2),chsel(3),:)),'Color',stcmap(2,:),'LineWidth',2);

% a(legsel(2)).LineStyle = '--';
% legend(a(legsel),legn)
xlim([4 38])
xlabel('frequency (Hz)')
ylabel('amplitude (a.u.)')
title('Simulated STN Spectra')
box off
grid on