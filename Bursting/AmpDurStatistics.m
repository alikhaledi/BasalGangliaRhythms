function BB = AmpDurStatistics(R,BB,condsel,CON)
BB = computeBetaBurstAmpDurStats_v2(R,BB);

ip(:,1) = [1 2 3]; % panel indice for plotting
%     else
%         condsel = [15 17 20];
%         ip(:,2) = [2 4 6];
%     end
%     BB = compute_BurstThreshold(R,BB,condsel,0);
%     R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STR->GPe','150% STR->GPe'};

subplot(2,3,1);
[h,l] = plotBurstAmplitudeHistogram(R,BB,condsel);
ylim(BB.plot.lims.burfreq(1,:))
p = get(gca,'Children');
%     if CON == 1
%         set(gca,'Children',p([1 2 3]));
%     else
%         set(gca,'Children',p([3 2 1]));
%     end
subplot(2,3,2);
plotBurstDurationHistogram(R,BB,condsel);
ylim(BB.plot.lims.burfreq(2,:))
p = get(gca,'Children');
p(1).FaceAlpha = 0.85;
%     if CON == 1
%         set(gca,'Children',p([1 2 3]));
%     else
%         set(gca,'Children',p([3 2 1]));
%     end
% subplot(2,3,3);
% plotBurstAmpDurScatter(R,BB,condsel)
%     p = get(gca,'Children');
%     set(gca,'Children',p([1 3 5 2 4 6]));
%     if CON == 1
%         set(gca,'Children',p([1 2 3]));
%     else
%         set(gca,'Children',p([3 2 1]));
%     end

%     cmap = brewermap(30,'Reds');
%     cntlevs = 0.5:0.5:15;
%     [h,l] = plotBurstAmplitudeDurationHistogram(R,BB,condsel,cmap,cntlevs);
%     cb = colorbar('southoutside');
%     cb.Position = [0.1548    0.0360    0.7589    0.0112];
%     cb.Label.String = 'Burst Rate (min^-1)';
%



featlist = {'ssAmp','ssDur','ssPPC','ssPow'};
featname = {'Burst Amplitude','Burst Duration','Burst CTX/STN Phase Sync.','Burst Power (Amp x Dur)'};
yl_list = {[0 180],[0 400],[0.5 1]}
for feat = 1:4
    subplot(2,3,feat+2)
    boundedline(linspace(10,190,19),BB.condstat.(featlist{feat})(1,:),BB.condstat.(featlist{feat})(4,:))
    
    xlabel('Connection Strength')
    ylabel(featname{feat})
    
    
   yyaxis right
   plot(linspace(10,190,19),[BB.segRate{:}],'k')
   ylabel('Burst Rate (min^-1)')
%     xlim([-1 0.75])
%     ylim(yl_list{feat})
    grid on
end




set(gcf,'Position',[313    45   1300    675])
%     figure
%     cmap = brewermap(30,'Blues');
%     cntlevs = 0.5:0.5:10;
%     [h,l] = plotBurstAmplitudeDurationHistogram(R,BB,[15 17 20],cmap,cntlevs);
%     cb = colorbar('southoutside');
%     cb.Position = [0.1548    0.0360    0.7589    0.0112];
%     cb.Label.String = 'Burst Rate (min^-1)';
%     set(gcf,'Position',[680    29   336   949])

