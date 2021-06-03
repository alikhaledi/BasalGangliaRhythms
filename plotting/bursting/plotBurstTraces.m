function plotBurstTraces(twin,m2Env,stnEnv,sw_twin,sw_PLV,befdPhi,middPhi,aftdPhi,dPhi,cmap)
% cmap = linspecer(3);
% cmap = brewermap(4,'RdYlBu');
oftab = [2e-7 10e-7 0.6 1; 2.5e-7 8e-7 0.6 1];
ostab = [1e-8 5e-8 0.05 0.01; 1e-8 2e-8 0.05 0.01];
lslist = {'-','--'};
npan = 5;
pvlist = nan(3,2,3);
for band = 1:2
    for state = 1:3
        subplot(2,npan,1+((band-1)*npan))
        normvar = mean(m2Env{band,state}(1,:)); %mean(m2Env{band,1}(:));
        %          normvar = mean(m2Env{band,1},2); % If plotting difference from baseline
        fy = mean(m2Env{band,state}-normvar,2);
        fz = std(m2Env{band,state}-normvar,[],2)./sqrt(size(m2Env{band,state},2));
        boundedline(twin{band,state},fy,fz,lslist{band},'cmap',cmap(state,:),'alpha','transparency',0.5); hold on
        if state>1
            base = (m2Env{band,4}-mean(m2Env{band,4}(1,:)));
            fy = (m2Env{band,state}-normvar);
            permSigStat(twin{band,state},base,fy,500,oftab(band,1),state*ostab(band,1),cmap(state,:),0.05,40)
        end
        title('M2 Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); % ylim([0 3]);
        grid on; box off; axis square
        
        subplot(2,npan,2+((band-1)*npan))
        normvar = mean(stnEnv{band,state}(1,:)); %mean(stnEnv{band,1}(:));
        fy = mean(stnEnv{band,state}-normvar,2);
        fz = std(stnEnv{band,state}-normvar,[],2)./sqrt(size(stnEnv{band,state},2));
        boundedline(twin{band,state},fy,fz,lslist{band},'cmap',cmap(state,:),'alpha','transparency',0.5); hold on
        if state>1
            base = (stnEnv{band,4} - mean(stnEnv{band,4}(1,:)));
            fy = (stnEnv{band,state}-normvar);
            permSigStat(twin{band,state},base,fy,500,oftab(band,2),state*ostab(band,2),cmap(state,:),0.05,40)
        end
        title('STN Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); %ylim([-1 3]);
        grid on; box off; axis square
        
        subplot(2,npan,3+((band-1)*npan))
        normvar = circ_mean(dPhi{band,state}(1,:),[],2);
        fy = circ_mean(dPhi{band,state},[],2)-normvar;
        fz = circ_std(dPhi{band,state},[],[],2)./sqrt(size(dPhi{band,state},2));
        boundedline(twin{band,state},wrapToPi(fy),fz,lslist{band},'cmap',cmap(state,:),'alpha','transparency',0.5); hold on
        if state>1
            base = wrapToPi(dPhi{band,4}-circ_mean(dPhi{band,4}(1,:),[],2));
            fy = wrapToPi(dPhi{band,state}-normvar);
            permSigStat(twin{band,state},base,fy,500,oftab(band,3),state*ostab(band,3),cmap(state,:),0.05,40)
        end
        title('M2/STN Phase'); ylabel('Relative Phase Angle'); xlabel('Time to Burst Onset (ms)'); %ylim([0.6 1.4]);
        grid on; box off; axis square
        
        subplot(2,npan,4+((band-1)*npan))
        normvar = mean(abs(sw_PLV{band,state}(1,:))); %mean(abs(sw_PLV{band,1}(:)));
        fy = mean(abs(sw_PLV{band,state}),2); %-normvar;
        fz = std(abs(sw_PLV{band,state}),[],2)./sqrt(size(sw_PLV{band,state},2));
        boundedline(sw_twin{band,state},fy,fz,lslist{band},'cmap',cmap(state,:),'alpha','transparency',0.5); hold on
        if state>1
            base = abs(sw_PLV{band,4});
            fy = abs(sw_PLV{band,state});
            permSigStat(sw_twin{band,state},base,fy,500,oftab(band,4),state*ostab(band,4),cmap(state,:),0.05,40)
        end
        title('M2/STN PLV'); ylabel('Normalized PLV'); xlabel('Time to Burst Onset (ms)'); ylim([0.6 1.1]);
        grid on; box off; axis square
        
        
        %         subplot(2,6,5+((band-1)*6))
        % %         p = polarhistogram(circ_mean(befdPhi{band,state},[],1),-pi:pi/32:pi,'Normalization','pdf'); hold on; p.FaceAlpha = 0.6; p.EdgeAlpha = 0;
        %         p(state) = polarhistogram(circ_mean(middPhi{band,state},[],1),-pi:pi/32:pi,'Normalization','probability');hold on;
        %         p(state).FaceAlpha = 0.6; p(state).EdgeAlpha = 0; p(state).FaceColor = cmap(state,:);
        %         title('Within Burst STN/M2 Phase')
        
        subplot(2,npan,5+((band-1)*npan))
        magvec = 0.1 + ((state-1).*0.025);
        
        % Radar plot
        Y = middPhi{band,state};
        X = middPhi{band,4};
        PS =circ_mean(Y,[],1);
        p(state) = plotPolarStem(PS,cmap(state,:)*1,magvec(1),'-','o',3,100);
        if state>1
            stattab(:,band,state) = radarStat(X,Y,cmap(state,:)*1,magvec(1),pi,'*',0.05);
        end

        title('Baseline Phase')
        rlim([0 0.16])
    end
end
l = legend(p,{'Fitted','A','B'});
l.Position = [0.9005 0.0582 0.0463 0.0661];
        a = gca;
a.FontSize = 16;
a.FontName = 'Arial';
function stattab = radarStat(X,Y,cmap,magvec,shift,symb,alpha)
[pv tabww] = circ_wwtest(circ_mean(X,[],1),circ_mean(Y,[],1));
stattab = [circ_mean(X(:))-circ_mean(Y(:)) circ_std(X(:))-circ_std(Y(:)) tabww{4,2} tabww{2,5} tabww{2,6}];
if pv<alpha
    polarscatter(circ_mean(Y(:)',[],2)+shift,magvec,symb,'MarkerEdgeColor',cmap)
end


%%% SCRIPT GRAVE
% % % Plot Within burst phase angles in radar plot + stat
% %         for i = 1:3
% %             switch i
% %                 case 1
% %                     X = befdPhi{band,1}; Y = befdPhi{band,state};
% %                 case 2
% %                     X = middPhi{band,1}; X2 = befdPhi{band,state}; Y = middPhi{band,state};
% %                 case 3
% %                     X = aftdPhi{band,1}; X2 = middPhi{band,state}; Y = aftdPhi{band,state};
% %             end
% %             PS =circ_mean(Y,[],1);
% %             p(i) = plotPolarStem(PS,cmap(state,:)*1,magvec(i),'-','o');
% %             if state>1
% %                 stattab(:,band,state) = radarStat(X,Y,cmap(state,:)*1,magvec(i),pi,'*',0.05);
% %                 if i>1
% %                     stattab2(:,band,state) = radarStat(X2,Y,cmap(state,:)*1,magvec(i),pi+0.2,'sqr',0.05);
% %                 end
% %             end
% %         end
% %         a = gca;
% %         a.RTick = [0.05 0.1 0.15];
% %         a.RTickLabel = {'Pre','Mid','Post'};
% %         title('Phase Progression')

