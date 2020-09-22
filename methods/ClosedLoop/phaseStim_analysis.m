%function phaseStim_analysis

close all
load tmpR
load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(1) '_stimM2_sensSTN_xsims.mat'],'xsimStore','uexsStore')
load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(1) '_stimM2_sensSTN_thresholdFitted.mat'],'intpow'); %,'burRP','plvStore'); %,'trajStore')
for band = 1:2
    
    site = 2;  baseCon = 1; %band = 2;
    [a maxind] = max(squeeze(intpow(site,band,2,:,baseCon)));
    [a minind] = min(squeeze(intpow(site,band,2,:,baseCon)));
    
%     pht = [maxind minind];
%     pht = [1 3 7];
        pht = [12 6];

    titname = {'Amplifying','Suppressive'};
    
    cm = {'r','b'};
    
    for phtype = 1:numel(pht)
        
        for stm = 1:2
            
            xsim_ip = xsimStore{stm,pht(phtype)};
            uexs = uexsStore{2,pht(phtype)};
            
            pU = uexs(R.IntP.phaseStim.sensStm(2),round(R.obs.brn*(1/R.IntP.dt))+1:end);
            indsS = SplitVec(find(abs(pU)>0),'consecutive'); % Split up data based upon the target threshold
            indsN = SplitVec(find(abs(pU)==0),'consecutive'); % Split up data based upon the target threshold
            
            segL = cellfun('length',indsS);
            % Now crop using minimum
            fsamp = (1/R.IntP.dt);
            burstSelection1 = find(segL>0); %((1.5/18)*fsamp));
            
            % Make sure the ends dont go off
            origin = cellfun(@(x) (x(1)-(0.3*fsamp)),indsS);
            burstSelection2 = find(origin>1);
            
            origin = cellfun(@(x) (x(end)+(1*fsamp)),indsS);
            burstSelection3 = find(origin<numel(pU));
            
            % Merge Selection
            burstSelection = intersect(burstSelection1,intersect(burstSelection2,burstSelection3));
            burstSelection(1) = []; % Remove first due to filter artefact
            %% Get the data features (bandpass
            X = xsim_ip([1 4],:);
            %     XF = ft_preproc_bandpassfilter(X,fsamp,[14 21],[],'fir');
            if band == 1
                fhz = 18; butf = [14 21];
            elseif band == 2
                fhz = 24; butf = [21 30];
            end
            BandWidth=3;
            SR =fsamp;
            wo = fhz/(SR/2);
            bw = BandWidth/(SR/2);
            qf =3;  % Q factor ( -qf);
            [bnt,ant] = iirnotch(wo,bw,qf);
%             [bpk,apk] = iirpeak(wo,bw,qf);
%                         [dum,bpk,apk] = ft_preproc_bandpassfilter(X, fsamp, butf,4,'but','twopass');
%                         XF =filtfilt(bpk,apk,X')';
            %
            % Alek's subtracted notch filter
            XF_nt =filtfilt(bnt,ant,X')';
            XF = X-XF_nt;
            
            for ch = 1:2
                XEnv(ch,:) = abs(hilbert(XF(ch,:)));
                XPhi(ch,:) = angle(hilbert(XF(ch,:)));
            end
            winsize(1) = 0.3*fsamp;
            winsize(2) = 1*fsamp;
            
            cycAmp{phtype,band,stm} = nan(numel(burstSelection),7);
            cycPhi{phtype,band,stm} = nan(numel(burstSelection),7);
            
            for seg = 1:numel(burstSelection)
                zeroP = indsS{burstSelection(seg)}(1);
                befWin = zeroP-winsize(1):zeroP;
                aftWin = zeroP:zeroP+winsize(2);
                
                befWav(:,seg) = XF(site,befWin);
                aftWav(:,seg) = XF(site,aftWin);
                
                befEnv{phtype,band,stm}(:,seg) = XEnv(site,befWin);
                aftEnv{phtype,band,stm}(:,seg) = XEnv(site,aftWin);
                
                befPhi{phtype,band,stm}(:,seg) = diff(unwrap(XPhi(:,befWin)))';
                %             befPhi(:,seg) = befPhi(:,seg)-befPhi(1,seg);
                aftPhi{phtype,band,stm}(:,seg) = diff(unwrap(XPhi(:,aftWin)))';
                %             aftPhi(:,seg) = aftPhi(:,seg)-aftPhi(1,seg);
                
                SWPLV{phtype,band,stm}(:,seg) = swPLV([befPhi{phtype,band,stm}(:,seg); aftPhi{phtype,band,stm}(:,seg)],0.15*fsamp);
                
                befEPhi(seg) = circ_mean(diff(unwrap(XPhi(:,befWin)))');
                aftEPhi(seg) = circ_mean(diff(unwrap(XPhi(:,aftWin)))');
                
                befdPhi(:,seg) = abs(diff(unwrap(XPhi(2,befWin))));
                aftdPhi(:,seg) = abs(diff(unwrap(XPhi(2,aftWin))));
                
                
                % Cycle by cycle
                %             X = pU(aftWin);
                X = XF(2,aftWin);
                zci = find(diff(sign(X))>1,6,'first'); % location of first 5 zero crossings
                if numel(zci)<6
                    zci(6-(6-numel(zci))+1:6) = nan;
                    warning('zci wrong')
                end
                zci = fix([zci(1)-(fsamp/18) zci zci(1)+(1*fsamp)]);
                zci = zeroP + zci;
                
                for cycle = 1:numel(zci)-1
                    winInds = zci(cycle):zci(cycle+1);
                    
                    Amp = XEnv(2,winInds);
                    Phi = diff(unwrap(XPhi(:,winInds)))';
                    %                 if cycle==1
                    cycAmp{phtype,band,stm}(seg,cycle) = mean(Amp);
                    cycPhi{phtype,band,stm}(seg,cycle) = circ_mean(Phi);
                    
                    %                 else
                    %                     cycAmp(seg,cycle) = mean(Amp-cycAmp(seg,1));
                    %                     cycPhi(seg,cycle) = circ_mean(Phi - cycPhi(seg,1));
                    %                 end
                end
            end
            
            % Across Trial Locking
            befPLV{phtype,band,stm} = splitapply(@PLV,befPhi{phtype,band,stm}',1:size(befPhi{phtype,band,stm},1));
            aftPLV{phtype,band,stm} = splitapply(@PLV,aftPhi{phtype,band,stm}',1:size(aftPhi{phtype,band,stm},1));
            
            % Slip Probability
            befPhSl = sum(befdPhi>prctile(befdPhi,75),2)/size(befdPhi,2);
            aftPhSl = sum(aftdPhi>prctile(aftdPhi,75),2)/size(aftdPhi,2);
        end
    end
end
a = 1;

cmap = brewermap(12,'Set1');
befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
TW = linspace(-winsize(1)/fsamp,winsize(2)/fsamp,size(SWPLV{band,stm},1));

% pht = [7 1];
%     pht = [1 5 7];

% First Plot the Amplitude Envelopes

for band = 1:2
    % mean Envelope
    subplot(3,2,sub2ind([2 3],band,1))    
    plot(befTW,mean(befEnv{phtype,band,1},2),'color','k','LineWidth',2);
        hold on
    plot(aftTW,mean(aftEnv{phtype,band,1},2),'color','k','LineWidth',2)
    
    for phtype = 1:numel(pht)
        plot(befTW,mean(befEnv{phtype,band,2},2),'color',cmap(pht(phtype),:),'LineWidth',2);
        plot(aftTW,mean(aftEnv{phtype,band,2},2),'color',cmap(pht(phtype),:),'LineWidth',2);
    end
    ylabel('STN Envelope'); xlabel('Time to stim onset (ms)')
    grid on; box off
    
    % Across Trial PLV
    subplot(3,2,sub2ind([2 3],band,2))    
    plot(befTW,befPLV{phtype,band,1},'color','k','LineWidth',2);
        hold on
    plot(aftTW,aftPLV{phtype,band,1},'color','k','LineWidth',2)
    
    for phtype = 1:numel(pht)
        plot(befTW,befPLV{phtype,band,2},'color',cmap(pht(phtype),:),'LineWidth',2);
        plot(aftTW,aftPLV{phtype,band,2},'color',cmap(pht(phtype),:),'LineWidth',2);
    end
    ylabel('Across-Trial PLV'); xlabel('Time to stim onset (ms)')
    grid on; box off
    
    
    % Within-Trial PLV
    subplot(3,2,sub2ind([2 3],band,3))    
    plot(TW,mean(SWPLV{phtype,band,1},2),'color','k','LineWidth',2);
        hold on
    plot(TW,mean(SWPLV{phtype,band,1},2),'color','k','LineWidth',2)
    
    for phtype = 1:numel(pht)
        plot(TW,mean(SWPLV{phtype,band,2},2),'color',cmap(pht(phtype),:),'LineWidth',2);
        plot(TW,mean(SWPLV{phtype,band,2},2),'color',cmap(pht(phtype),:),'LineWidth',2);
    end    
    ylabel('Within-Trial PLV'); xlabel('Time to stim onset (ms)')
    grid on; box off
    
    
    
end




function coeff = PLV(phi)
coeff = abs(sum(exp(1i*phi)))./numel(phi);
end
function out = swPLV(phi,winsize)
slidePhi = slideWindow(phi,winsize,winsize*0.95);
out = splitapply(@PLV,slidePhi,1:size(slidePhi,2));
end
function y=rodd(x)
y = 2*floor(x/2)+1;
end