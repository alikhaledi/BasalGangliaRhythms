clear; close all
load tmpR
load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(1) '_stimM2_sensSTN_xsims.mat'],'xsimStore','uexsStore')
pht = [1 6];
titname = {'Amplifying','Suppressive'};
for phtype = 1:2
    
    for base = 2; %1:2
        if base == 1
            cm = 'r';
        elseif base == 2
            cm = 'b';
        end
        xsim_ip = xsimStore{base,pht(phtype)};
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
        
        origin = cellfun(@(x) (x(end)+(0.6*fsamp)),indsS);
        burstSelection3 = find(origin<numel(pU));
        
        % Merge Selection
        burstSelection = intersect(burstSelection1,intersect(burstSelection2,burstSelection3));
        burstSelection(1) = []; % Remove first due to filter artefact
        %% Get the data features (bandpass
        X = xsim_ip([1 4],:);
        %     XF = ft_preproc_bandpassfilter(X,fsamp,[14 21],[],'fir');
        fhz = 18;
        BandWidth=3;
        SR =fsamp;
        wo = fhz/(SR/2);
        bw = BandWidth/(SR/2);
        qf =3;  % Q factor ( -qf);
        [bpk,apk] = iirnotch(wo,bw,qf);
        %         [bpk,apk] = iirpeak(wo,bw,qf);
        XF =filtfilt(bpk,apk,X')';
        XF = X-XF;
        figure(100)
        
        for ch = 1:2
            XEnv(ch,:) = abs(hilbert(XF(ch,:)));
            XPhi(ch,:) = angle(hilbert(XF(ch,:)));
        end
        winsize(1) = 0.3*fsamp;
        winsize(2) = 0.6*fsamp;
        
        cycAmp = nan(numel(burstSelection),6);
        cycPhi = nan(numel(burstSelection),6);
        for seg = 1:numel(burstSelection)
            zeroP = indsS{burstSelection(seg)}(1);
            befWin = zeroP-winsize(1):zeroP;
            aftWin = zeroP:zeroP+winsize(2);
            
            X = pU(aftWin);
            
             zci = find(diff(sign(X))>1,6,'first'); % location of first 5 zero crossings
            zci = fix([zci(1)-(fsamp/18) zci]);
            zci = zeroP + zci;
            
            for cycle = 1:numel(zci)-1
                winInds = zci(cycle):zci(cycle+1);
            
                cycAmp(seg,cycle) = mean(XEnv(2,winInds));
                cycPhi(seg,cycle) = circ_mean(diff(unwrap(XPhi(:,winInds)))');
                
                if cycle>1
                    cycAmp(seg,cycle) = cycAmp(seg,cycle)-cycAmp(seg,1);
                 cycPhi(seg,cycle) =  cycPhi(seg,cycle)- cycPhi(seg,1);
                end
            end
        end
        
        %     rose(aftPhi,32);
                % Slip Probability
        befPhSl = sum(befdPhi>prctile(befdPhi,75),2)/size(befdPhi,2);
        aftPhSl = sum(aftdPhi>prctile(aftdPhi,75),2)/size(aftdPhi,2);
        
        subplot(3,2,phtype)
        befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
        plot(befTW,mean(befEnv,2),'color',cm);
        hold on
        %     plot(befTW,mean(befWav,2));
        
        aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
        plot(aftTW,mean(aftEnv,2),'color',cm)
        %         plot(aftTW,mean(aftWav,2))
        
        xlabel('Time to Stim Onset (ms)'); ylabel('Beta Envelope')
        %     ylim([0.8 1.3]*1e-7)
        title(titname{phtype})
        
        subplot(3,2,phtype+2)
        befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
        plot(befTW(1:end),wrapTo2Pi(circ_mean(befPhi,[],2)),'color',cm);
        hold on
        aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
        plot(aftTW(1:end),wrapTo2Pi(circ_mean(aftPhi,[],2)),'color',cm)
        %     ylim([-0.5 0.75])
        xlabel('Time to Stim Onset (ms)'); ylabel('CTX-STN Phase')
        
        % Phase Distribution
        %     subplot(2,2,phtype+2)
        %     rose(befPhi,32);yeah
        %     hold on

        
        subplot(3,2,phtype+4)
        befTW = linspace(-winsize(1)/fsamp,0,winsize(1)+1);
        plot(befTW(2:end),befPhSl,'color',cm);
        hold on
        aftTW = linspace(0,winsize(2)/fsamp,winsize(2)+1);
        plot(aftTW(2:end),aftPhSl,'color',cm)
        xlabel('Time to Stim Onset (ms)'); ylabel('Phase Slip Probability')
        %     ylim([0.05 0.175])
        
        
        figure(200)
        subplot(1,2,phtype)
        rose(befEPhi); hold on; rose(aftEPhi)
        title(titname{phtype})
    end
end
[pval table] = circ_wwtest(befEPhi,aftEPhi)
