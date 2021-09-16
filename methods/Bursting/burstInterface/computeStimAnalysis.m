function [R,BB] = computeStimAnalysis(R,fresh)

rootan = [R.rootn 'data\phaseLockedStim'];
close all;
conname = {'HD','','PS'};
statename = {'Fitted','A','B'};
SScomb = 1;
CON = 1;

if fresh
    % Load in Data
    BB = [];
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],'feat_sim_save');
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip');
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_Rout' num2str(SScomb) '.mat'],'Rout');
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_pU_save' num2str(SScomb) '.mat'],'pU_save');
    fsamp = 1/Rout.IntP.dt;
    for stm = 1:2
        % Set the phase range
        if stm==1
            philist = 1;
        elseif stm ==2
            philist = 1:12;
        end
        % Do base model
        state = 1;
        Xbase = xsim_ip{1,state}{1};
        Xpair = [1 4]; % The pair of pop indices to be analysed
        for phi = philist
            % Get Data
            X = xsim_ip{stm,state}{phi}{1};
            pU{phi,stm} = pU_save{state}{phi};
            
            % Get Burst Inds
            stimSelInds = SplitVec(find(abs(pU{phi,stm})>0),'consecutive'); % Split up data based upon the target threshold
            
            % Setup window in which to look
            winsize(1) = 0.075*fsamp;
            winsize(2) = 0.075*fsamp; %0.15*fsamp;
            winsize(3) = 0.075*fsamp;
            winsize(4) = 0.075*fsamp;
            winloc(1) =  -0.35*fsamp;
            winloc(2) =  0*fsamp;
            winloc(3) =  +0.5*fsamp;
            
            stimSelInds = cropBurstSelection(stimSelInds,winsize,size(X,2));% Remove bursts at ends
            
            XS{phi,1} = Xbase{1}(Xpair,[stimSelInds{:}]); % No stim
            if stm>1
                XS{phi,2} = X(Xpair,[stimSelInds{:}]); % with stim
            end
            
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
                
                % Get Filtered/Analytic Signal
                [~,XF,XEnv,XPhi,epsAmp,segL{band,phi,stm},segA{band,phi,stm}] = simpleBurstDefine(X,fsamp,frqz,minper); % Remember no burst definition done here, just used for filter
                
                % Get burst overlap
                [statBurstOverl{band,phi,stm},sampBurstOverl{band,phi,stm}] = burstOverlap(XF,fsamp,frqz,minper,band);
                
                % Compute Connectivity Matrix
                [Pw pHz] = pwelch(X',fsamp,[],fsamp,fsamp);
                specMat{band,phi,stm} = Pw(pHz>4 & pHz<=48,:);
                [connectMat{band,phi,stm},diffConnectCI{band,phi,stm}] = computePLVConMat(XPhi(:,[stimSelInds{:}]),band);
                
                [twin{band,phi,stm},m2Env{band,phi,stm},stnEnv{band,phi,stm},...
                    dPhi{band,phi,stm},sw_twin{band,phi,stm},sw_PLV{band,phi,stm},...
                    aftdPhi{band,phi,stm},befdPhi{band,phi,stm},middPhii{band,phi,stm},allEnv{band,phi,stm}] = getBurstTraces(stimSelInds,winsize,winloc,fsamp,XEnv,XPhi,Xpair);
                
                % Compute relative Burst Onsets
                initTime = relativePeakTiming(allEnv{band,phi,stm},twin{band,phi,stm},epsAmp,minper,frqz,fsamp);
            end
        end
    end
    save([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysis_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
        'connectMat','specMat','statBurstOverl','twin','m2Env','stnEnv','dPhi','sw_twin','sw_PLV',...
        'aftdPhi','befdPhi','middPhii','Rout','XS','allEnv','pU','segL','segA','diffConnectCI')
else
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysis_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
        'connectMat','specMat','statBurstOverl','twin','m2Env','stnEnv','dPhi','sw_twin','sw_PLV',...
        'aftdPhi','befdPhi','middPhii','Rout','XS','allEnv','pU','segL','segA','diffConnectCI')
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip');
end


Hz = Rout.frqz;
for stm = 1:2
    for phi = 1:12
        if stm == 1
            phiEf = 1;
        else
            phiEf = phi;
        end
        %         spec = spectraSave{phiEf,stim};
        [F,Hz] = pwelch(XS{phiEf,stm}',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
        [C,Hz] = mscohere(XS{phiEf,stm}(1,:)',XS{phiEf,stm}(2,:)',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
        
        spec = [F C];
        powspec_save(:,:,stm,phi,1) = spec;
        intpow(:,1,stm,phi,1) = sum(spec(Hz>14 & Hz<=21,:))./numel(find(Hz>14 & Hz<=21));
        maxpow(:,1,stm,phi,1) = max(spec(Hz>14 & Hz<=21,:));
        intpow(:,2,stm,phi,1) = sum(spec(Hz>21 & Hz<=30,:))./numel(find(Hz>21 & Hz<=30));
        maxpow(:,2,stm,phi,1) = max(spec(Hz>21 & Hz<=30,:));
        intpow(:,3,stm,phi,1) = sum(spec(Hz>14 & Hz<=30,:))./numel(find(Hz>14 & Hz<=30));
        maxpow(:,3,stm,phi,1) = max(spec(Hz>14 & Hz<=30,:));
    end
end

figure
cmap = brewermap(12,'Set1');
phangle = [6 12];
Rout.cmap = [cmap(phangle(1)-1,:); cmap(phangle(2)-1,:); 0 0 0];
plotStimExampleTrace(pU,xsim_ip,phangle,Rout)
set(gcf,'Position',[747   316   813   489])

phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12

figure; closedLoopSpectralPlot(R,phaseShift,1,intpow,powspec_save,Hz)
figure; plotStimTraces2(twin,m2Env,stnEnv,sw_twin,sw_PLV,dPhi,befdPhi,middPhii,aftdPhi)
cmap = brewermap(12,'Set1');
philist = [100 6 12];
cmap = [0 0 0; cmap(philist(2)-1,:); cmap(philist(3)-1,:)];
set(gcf,'Position',[401         321        1368         657])

figure; bplotStimBurstLength(segA,segL,'sda',cmap,2)
% end