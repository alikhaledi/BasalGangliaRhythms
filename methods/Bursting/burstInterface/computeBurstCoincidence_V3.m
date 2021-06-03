function [R,BB] = computeBurstCoincidence_V3(R)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
close all;
conname = {'HD','','PS'};
statename = {'Fitted','A','B'};

for CON = [1 3]
    BB = [];
    load([rootan '\BB_' R.out.tag '_DiscreteData.mat'],'dataSelect','dataProperties')
    
    for state = 1:3
        %% Find Burst Epochs
        % Get Data
        X = dataSelect{CON}{state}{1};
        for band = 1:2
            if band == 1
                fhz = 18; butf = [14 21];
            elseif band == 2
                fhz = 24; butf = [21 30];
            end
            
            % Set parameters
            fsamp = 1/R.IntP.dt;
            frqz = butf;
            minper = 3; % Min 3 cycles
            %         peakHz = dataProperties(4,1,CON); % This loads in the peak STN frequency
            [burstSelInds,XF,XEnv,XPhi] = simpleBurstDefine(X,fsamp,frqz,minper);
            
            
            
            % Setup window in which to look
            winsize(1) = 0.3*fsamp;
            winsize(2) = 1*fsamp;
            burstSelInds = cropBurstSelection(burstSelInds,winsize,size(X,2));% Remove bursts at ends
            
            [twin{band,state},m2Env{band,state},stnEnv{band,state},...
             dPhi{band,state},sw_twin{band,state},sw_PLV{band,state},...
             aftdPhi{band,state},aftdPhi2{band,state}] = getBurstTraces(burstSelInds,winsize,fsamp,XEnv,XPhi);
            
        end
    end
    for band = 1:2
        for state = 1:3
            figure(CON)
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
    figure(CON)
    set(gcf,'Position',[680         367        1036         611])
    l = legend({'Fitted','A','B'});
    l.Position = [0.9012 0.0360 0.0753 0.0769];
end

a = 1;
% end