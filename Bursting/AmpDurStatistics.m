function BB = AmpDurStatistics(R,BB,condsel,CON)
BB = computeBetaBurstAmpDurStats_v2(R,BB);

ip(:,1) = [1 2 3]; % panel indice for plotting
%     else
%         condsel = [15 17 20];
%         ip(:,2) = [2 4 6];
%     end
%     BB = compute_BurstThreshold(R,BB,condsel,0);
%     R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STR->GPe','150% STR->GPe'};

subplot(3,2,1);
[h,l] = plotBurstAmplitudeHistogram(R,BB,condsel);

subplot(3,2,2);
plotBurstDurationHistogram(R,BB,condsel);
% subplot(2,3,3);
% plotBurstAmpDurScatter(R,BB,condsel)

featlist = {'ssAmp','ssDur','ssPPC','ssPow'};
featname = {'Burst Amplitude','Burst Duration','Burst CTX/STN Phase Sync.','Burst Power (Amp x Dur)'};
yl_list = {[0 20],[0 250],[0 0.5],[0 2500]}
for feat = 1:4
    subplot(3,2,feat+2)
    boundedline(linspace(10,190,19),BB.condstat.(featlist{feat})(1,:),BB.condstat.(featlist{feat})(4,:))
    
    xlabel('Connection Strength')
    ylabel(featname{feat})
    
    
   yyaxis right
   plot(linspace(10,190,19),[BB.segRate{:}],'k')
   ylabel('Burst Rate (min^-1)')
%     xlim([-1 0.75])
   yyaxis left

    ylim(yl_list{feat})
    grid on
end



set(gcf,'Position',[952    30   840   948])
% set(gcf,'Position',[313    45   1300    675])
%     figure
%     cmap = brewermap(30,'Blues');
%     cntlevs = 0.5:0.5:10;
%     [h,l] = plotBurstAmplitudeDurationHistogram(R,BB,[15 17 20],cmap,cntlevs);
%     cb = colorbar('southoutside');
%     cb.Position = [0.1548    0.0360    0.7589    0.0112];
%     cb.Label.String = 'Burst Rate (min^-1)';
%     set(gcf,'Position',[680    29   336   949])

