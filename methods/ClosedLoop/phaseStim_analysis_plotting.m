
subplot(3,2,phtype)
befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
plot(befTW,mean(befEnv{band,stm},2),'color',cm{stm});
hold on
%     plot(befTW,mean(befWav,2));
aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
plot(aftTW,mean(aftEnv{band,stm},2),'color',cm{stm})
%         plot(aftTW,mean(aftWav,2))

xlabel('Time to Stim Onset (ms)'); ylabel('Beta Envelope')
ylim([3 25]*1e-7)
title(titname{phtype})

% Time series of phase angle
%             subplot(3,2,phtype+2)
%             befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
%             plot(befTW(1:end),wrapTo2Pi(circ_mean(befPhi{band,stm},[],2)),'color',cm{stm});
%             hold on
%             aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
%             plot(aftTW(1:end),wrapTo2Pi(circ_mean(aftPhi{band,stm},[],2)),'color',cm{stm})
%             %     ylim([-0.5 0.75])
%             xlabel('Time to Stim Onset (ms)'); ylabel('CTX-STN Phase')

subplot(3,2,phtype+2)
TW = linspace(-winsize(1)/fsamp,winsize(2)/fsamp,size(SWPLV{band,stm},1));
plot(TW,mean(SWPLV{band,stm},2),'color',cm{stm});
hold on
%     ylim([-0.5 0.75])
xlabel('Time to Stim Onset (ms)'); ylabel('Within Trial PLV')

subplot(3,2,phtype+4)
befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
plot(befTW(1:end),befPLV{band,stm},'color',cm{stm});
hold on
aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
plot(aftTW(1:end),aftPLV{band,stm},'color',cm{stm})
xlabel('Time to Stim Onset (ms)'); ylabel('Across Trial PLV')
ylim([0.4 1])

figure(300)
subplot(2,2,phtype)
A = cycAmp{band,stm}; %nanmean(cycAmp,1);
A = 100*(A-mean(A(:,1),1))/mean(A(:,1),1);
plot(1:7,mean(A))
errorbar(1:7,mean(A),std(A)./sqrt(numel(A)),'Color',cm{stm});
hold on
a = gca;
a.XTick = 1:7;
a.XTickLabel = {'-1 cycle','1st','2nd','3rd','4th','5th','+600ms'}
ylabel('% Change in Amplitude')
title(titname{phtype})

subplot(2,2,phtype+2)
P = cycPhi{band,stm};
P = P-circ_mean(P(:,1));
P = wrapToPi(P);
%     plot(1:7,P)
errorbar(1:7,circ_mean(P),circ_std(P)./sqrt(numel(P)),'Color',cm{stm});
hold on
a = gca;
a.XTick = 1:7;
a.XTickLabel = {'-1 cycle','1st','2nd','3rd','4th','5th','+600ms'};
ylabel('change in STN/M2 Phase')

figure(400)
subplot(1,2,phtype)
A = cycAmp{stm}; %nanmean(cycAmp,1);
p =  polarplot(circ_mean(P),1:7);
p.Marker = '>';
p.MarkerFaceColor = cm{stm};
p.Color = cm{stm};
thetalim([0-35 0+45])
hold on
title(titname{phtype})
ax = gca;
ax.RAxis.Label.String = 'Cycle Number';
ax.RAxis.TickValues = 1:7;
ax.RAxis.TickLabels = {'-1 cycle','1st','2nd','3rd','4th','5th','+600ms'};
ax.ThetaAxis.Label.String = 'STN/M2 Relative Phase';

% Find window in which there is a maximum effect
A1 = mean(aftEnv{band,1},2);
A2 = mean(aftEnv{band,2},2);
[dum maxind] = max(abs(A1-A2));
maxwin = [maxind-(0.15*fsamp):maxind+(0.15*fsamp)];
maxwin = maxwin(maxwin>0);
%         splitapply(@PLV,X,1:size(X,2))
for stm = 1:2
    figure(600+band)
    subplot(1,3,phtype)
    X = aftPhi{band,stm}(maxwin,:);
    p = polarhistogram(X,-pi:pi/16:pi,'Normalization','probability');
    p.FaceColor = cm{stm};
    hold on
    %             polarplot([0 sum(exp(i*X(:)))/numel(X(:))]*0.1)
    title(titname{phtype})
    plvlist(stm,phtype) = abs(sum(exp(1i*X(:))))/numel(X(:));
end

figure(600+band)
subplot(1,3,3)
b = bar(diff(plvlist));
a =gca;
a.XTickLabel = titname;
ylabel('Change in PLV ')
set(gcf,'Position',[543         542        1141         436])

% This will be for the Suppressive Phase Only
figure(800)
bcl = {'g','c'};
for band = 1:2
    subplot(2,1,1)
    befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
    A = mean(befEnv{band,2},2)-mean(befEnv{band,1},2);
    plot(befTW,A,'color',bcl{band});
    hold on
    %     plot(befTW,mean(befWav,2));
    aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
    A = mean(aftEnv{band,2},2)-mean(aftEnv{band,1},2);
    plot(aftTW,A,'color',bcl{band})
    xlabel('Time to Stim Onset (ms)')
    ylabel('Mean Envelope')
    
    subplot(2,1,2)
    befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
    A = befPLV{band,2}-befPLV{band,1};
    plot(befTW,A,'color',bcl{band});
    hold on
    %     plot(befTW,mean(befWav,2));
    aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
    A = aftPLV{band,2}-aftPLV{band,1};
    p(band) = plot(aftTW,A,'color',bcl{band})
    xlabel('Time to Stim Onset (ms)')
    ylabel('Within Trial PLV')
    legend(p,{'B1','B2'})
    
    
    subplot(2,1,3)
    TW = linspace(-winsize(1)/fsamp,winsize(2)/fsamp,size(SWPLV{band,stm},1));
    plot(TW,mean(SWPLV{band,stm},2),'color',cm{stm});
    hold on
    %     ylim([-0.5 0.75])
    xlabel('Time to Stim Onset (ms)'); ylabel('Within Trial PLV')
    
    
    
end
set(gcf,'Position',[924   205   673   773])
[pval table] = circ_wwtest(befEPhi,aftEPhi)
