function plotBurstTraces(twin,m2Env,stnEnv,sw_twin,sw_PLV,aftdPhi)
for band = 1:2
    for state = 1:3
        normvar = mean(m2Env{band,1}(:));
        subplot(2,4,1+((band-1)*4))
        plot(twin{band,state},mean(m2Env{band,state},2)./normvar); hold on;
        title('M2 Envelope'); ylim([0 3]); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)');
        grid on; box off;
        
        subplot(2,4,2+((band-1)*4))
        normvar = mean(stnEnv{band,1}(:));
        plot(twin{band,state},mean(stnEnv{band,state},2)./normvar); hold on;
        title('STN Envelope'); ylim([0 3]); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)');
        grid on; box off;
        
        subplot(2,4,3+((band-1)*4))
        normvar = mean(sw_PLV{band,1}(:));
        plot(sw_twin{band,state},mean(sw_PLV{band,state},2)./normvar); hold on;
        title('M2/STN PLV'); ylim([0.6 1.4]); ylabel('Normalized PLV'); xlabel('Time to Burst Onset (ms)');
        grid on; box off;
        
        subplot(2,4,4+((band-1)*4))
        p = polarhistogram(aftdPhi{band,state}(:),-pi:pi/64:pi,'Normalization','pdf'); hold on
        title('STN/M2 Relative Phase')
        p.FaceAlpha = 0.6;
        p.EdgeAlpha = 0;
    end
end
l = legend({'Fitted','A','B'});
l.Position = [0.9012 0.0360 0.0753 0.0769];
