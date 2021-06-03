function plotStimTraces(twin,m2Env,stnEnv,sw_twin,sw_PLV,dPhi,befdPhi,middPhi,aftdPhi)
cmap = brewermap(12,'Set1');
ls = {'-','-','-'};
% cmap = brewermap(4,'RdYlBu');
philist = [100 1 7];
cmap = [0 0 0; cmap(philist(2),:); cmap(philist(3),:)];

for band = 1:2
    for phi = 1:3
        if philist(phi) == 100
            stm = 1;
            pheff =  1;
        else
            stm = 2;
            pheff =  philist(phi);
        end
        figure(1)
        subplot(2,6,1+((band-1)*6))
        normvar = mean(m2Env{band,pheff,stm}(1,:)); %mean(m2Env{band,1}(:));
        fy = mean(m2Env{band,pheff,stm}-normvar,2);
        fz = std(m2Env{band,pheff,stm}-normvar,[],2)./sqrt(size(m2Env{band,pheff,stm},2));
        boundedline(twin{band,pheff,stm},fy,fz,'cmap',cmap(state,:)); hold on
%         plot(twin{band,state},fy,'Color',cmap(state,:)); hold on;
%         if state>1
%             base = (m2Env{band,1}-normvar);
%             fy = (m2Env{band,state}-normvar);
%             permSigStat(twin{band,state},base,fy,200,oftab(band,1),state*ostab(band,1),cmap(state,:))   
%         end
        title('M2 Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); % ylim([0 3]);
        grid on; box off;
        
        subplot(2,6,2+((band-1)*6))
        normvar = mean(stnEnv{band,pheff,stm}(1,:)); %mean(stnEnv{band,1}(:));
        fy = mean(stnEnv{band,pheff,stm}-normvar,2);
        fz = std(stnEnv{band,pheff,stm}-normvar,[],2)./sqrt(size(stnEnv{band,pheff,stm},2));
        boundedline(twin{band,pheff,stm},fy,fz,'cmap',cmap(state,:)); hold on
%         if state>1
%             base = (stnEnv{band,1}-normvar);
%             fy = (stnEnv{band,state}-normvar);
%             permSigStat(twin{band,state},base,fy,200,oftab(band,2),state*ostab(band,2),cmap(state,:))   
%         end        
        title('STN Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); %ylim([-1 3]); 
        grid on; box off;
        
        subplot(2,6,3+((band-1)*6))
        fy = circ_mean(dPhi{band,pheff,stm},2);
        fz = circ_std(dPhi{band,pheff,stm},[],2)./sqrt(size(dPhi{band,pheff,stm},2));
        boundedline(twin{band,pheff,stm},fy,fz,'cmap',cmap(state,:)); hold on
        title('STN Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); %ylim([-1 3]); 
        grid on; box off;
        
        subplot(2,6,4+((band-1)*6))
        normvar = mean(abs(sw_PLV{band,pheff,stm}(1,:))); %mean(abs(sw_PLV{band,1}(:)));
        fy = mean(abs(sw_PLV{band,pheff,stm}),2);
        fz = std(abs(sw_PLV{band,pheff,stm}),[],2)./sqrt(size(sw_PLV{band,pheff,stm},2));
        boundedline(sw_twin{band,pheff,stm},fy,fz,'cmap',cmap(state,:)); hold on
        title('M2/STN PLV'); ylabel('Normalized PLV'); xlabel('Time to Burst Onset (ms)'); %ylim([0.6 1.4]); 
        grid on; box off;
        
        subplot(2,6,5+((band-1)*6))
%         p = polarhistogram(circ_mean(befdPhi{band,state},[],1),-pi:pi/32:pi,'Normalization','pdf'); hold on; p.FaceAlpha = 0.6; p.EdgeAlpha = 0;
        p(state) = polarhistogram(circ_mean(middPhi{band,pheff,stm},[],1),-pi:pi/32:pi,'Normalization','probability');hold on;
        p(state).FaceAlpha = 0.6; p(state).EdgeAlpha = 0; p(state).FaceColor = cmap(state,:);
        title('Within Burst STN/M2 Phase')
        
        subplot(2,6,6+((band-1)*6))
        magvec = [0.05 0.1 0.15] + ((state-1).*0.01);
        X =circ_mean(befdPhi{band,pheff,stm},[],1);
        plotPolarStem(X,cmap(state,:)*1,magvec(1));
        X =circ_mean(middPhi{band,state},[],1);
        plotPolarStem(X,cmap(state,:)*1,magvec(2));
        X =circ_mean(aftdPhi{band,state},[],1);
        plotPolarStem(X,cmap(state,:)*1,magvec(3));
        
        title('Phase Progression')
    end
end
l = legend(p,{'Fitted','Suppressive','Promoting'});
l.Position = [0.9005 0.0582 0.0463 0.0661];
