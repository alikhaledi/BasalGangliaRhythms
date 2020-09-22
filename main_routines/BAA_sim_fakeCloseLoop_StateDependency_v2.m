function [R] = BAA_sim_fakeCloseLoop_StateDependency_v2(Rorg,modID,simtime,fresh)
% Load in model data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v3(R,modID,simtime);]
% or load from preload
load([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
R.filepathn =  Rorg.filepathn;
R.rootn =  Rorg.rootn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')
% pause(1)

%% Define stimulation conditions
% Stimualating M2
senssite = 4; % STN
stimsite = 1; % M2
stim_sens = 'stimM2_sensSTN';

% Stimulating  STN
% senssite = 1; % M2
% stimsite = 4; % STN
% stim_sens = 'stimSTN_sensM2';

R.IntP.phaseStim.sensStm = [senssite stimsite];

P = permMod{1}.par_rep{1};
% Simulation Coniditions
R.obs.csd.df = 0.5;
R = setSimTime(R,simtime);
R.obs.trans.norm = 0; % No normalization of spectra
R.obs.brn = R.obs.brn;
R.obs.gainmeth{1} = 'STIM';
R.IntP.intFx = @spm_fx_compile_120319_stim;
phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12
tplot = 0;
% Give all timeseries the same input - makes comparable
rng(5453)
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);

if fresh
    for ctype = 1
        %% CLEAR
        intpow = [];
        maxpow = [];
        
        burRate = [];
        burdur = [];
        burAmp = [];
        burAmpMid = [];
        burPPC = [];
        burInt = [];
        
        durStore = [];
        ampStore = [];
        ppcStore = [];
        siStore = [];
        
        rootan = [Rorg.rootn 'data\rat_InDirect_ModelComp\ConnectionSweep'];
        
        if ctype == 1 % M2 STN
            CON = 1;
        elseif ctype == 2 % GPe STN
            CON = 3;
        end
        
        load([rootan '\BB_'  R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
        
        NcS = ck_1(CON,1:2:end);
        NcS = [1 ck_1(CON,:)];
        
        for cond = 1:numel(NcS)
            
            % Get Base Parameters
            Pbase = P;
            
            % Modulate connectivity (state dependency)
            if ctype ==1
                Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*NcS(cond)); %
            elseif ctype == 2
                Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*NcS(cond)); %
            end
            %     intpow = []; maxpow = [];
            for p = 1:numel(phaseShift) %[1 10] %
                
                % Initialize variables
                uc_ip = {}; feat_sim = {}; xsim_ip = {};
                uc_ip{1} = uc;
                
                %% Setup stim parameterss
                switch stim_sens
                    case 'stimM2_sensSTN'
                R = typeIstimPars_v3(R);
                    case 'stimSTN_sensM2'
                R = typeIIstimPars(R);
                end
                % Modulate the phase
                R.IntP.phaseStim.phaseshift = phaseShift(p);
                
                %% Simulate Base Model
                R.IntP.phaseStim.switch = 0 ;
                [~,~,feat_sim{1},xsim_gl,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
                
                %% Work out the threshold
                [~,R] = zeroCrossingPhaseStim_v3([],R,[],xsim_gl{1},R.IntP.dt);
               eps(p) =  R.IntP.phaseStim.eps;
                %% Resimulate with Phase-Locked Input
                R.IntP.phaseStim.switch = 1;
                
                % Simulate with Stimulation
                [~,~,feat_sim{2},~,xsim_ip{2},~,Rout]  = computeSimData(R,m,uc_ip{1},Pbase,0);
                tplot = 0;
                    load([R.rootn 'data\rat_InDirect_ModelComp\phaseStimSave\stim_tmp'],'uexs')
                if tplot == 1
                    pU = uexs(4,round(R.obs.brn*(1/R.IntP.dt))+1:end);
                    % Optional Plots for TimeSeries
                    figure(100)
                    clf
                    a(1) = subplot(3,1,1);
                    plot(Rout.IntP.tvec_obs,pU);
                    
                    a(2) = subplot(3,1,2);
                    plot(Rout.IntP.tvec_obs,xsim_ip{1}{1}(1,2:end));
                    hold on
                    plot(Rout.IntP.tvec_obs,xsim_ip{2}{1}(1,2:end));
                    
                    a(3) = subplot(3,1,3);
                    plot(Rout.IntP.tvec_obs,xsim_ip{1}{1}(4,2:end)+5e-6);
                    hold on
                    plot(Rout.IntP.tvec_obs,xsim_ip{2}{1}(4,2:end)+5e-6);
                    linkaxes(a,'x')
                    %                     xlim([24.75 25.75])
                    R.plot.outFeatFx({feat_sim{1}},{feat_sim{2}},R.data.feat_xscale,R,1,[])
                end
                
                % Re-compute bursts (simulated data)
                % Low Beta
                R.condname = {'1','2'};
                R.BB.powfrq = 18;
                R.BB.cohfrq = 18;
                [R,BB] = compute_BetaBursts_Simulated(R,{xsim_ip{1}{1} xsim_ip{2}{1}});
                R.BB.thresh_prctile = 75;% o85; tl 80
                R.BB.threshold_type = 'baseModelThresh';
                BB = compute_BurstThreshold(R,BB,1,0);
                
                R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
                BB.plot.durlogflag = 0;
                R.BB.pairInd = [1 4]; % KEEP LOOKING AT STN BURSTS!
                BB = defineBetaEvents(R,BB);
                
                % High Beta
                R.BB.powfrq = 24;
                R.BB.cohfrq = 24;
                [R,BBH] = compute_BetaBursts_Simulated(R,{xsim_ip{1}{1} xsim_ip{2}{1}});
                BBH = compute_BurstThreshold(R,BBH,1,0);
                BBH = defineBetaEvents(R,BBH);
                %% Now get statistics of the simulations (spectra/burst features)
                
                for stm = 1:2
                    spec = [squeeze(feat_sim{stm}(1,1,1,1,:)) squeeze(feat_sim{stm}(1,4,4,1,:)) squeeze(feat_sim{stm}(1,4,1,4,:))];
                    powspec_save(:,:,stm,p,cond) = spec;
                    intpow(:,1,stm,p,cond) = sum(spec(R.frqz>14 & R.frqz<=21,:));
                    maxpow(:,1,stm,p,cond) = max(spec(R.frqz>14 & R.frqz<=21,:));
                    intpow(:,2,stm,p,cond) = sum(spec(R.frqz>21 & R.frqz<=30,:));
                    maxpow(:,2,stm,p,cond) = max(spec(R.frqz>21 & R.frqz<=30,:));
                    intpow(:,3,stm,p,cond) = sum(spec(R.frqz>14 & R.frqz<=30,:));
                    maxpow(:,3,stm,p,cond) = max(spec(R.frqz>14 & R.frqz<=30,:));
                    
                    xsimStore{stm,p} = xsim_ip{stm}{1};
                    uexsStore{stm,p} = uexs;
                    
                    
                    if numel(BB.segAmp{2})>2
                        burRate(:,stm,p,cond) = percentageChange(BB.segRate,stm);
                        burdur(:,stm,p,cond) = [nanmedian(percentageChange(BB.segDur,stm)) npCI(percentageChange(BB.segDur,stm))];
                        burAmp(:,stm,p,cond) = [nanmedian(percentageChange(BB.segAmp,stm)) npCI(percentageChange(BB.segAmp,stm))];
                        burAmpMid(:,stm,p,cond) = [nanmedian(percentageChange(BB.segAmpMid,stm)) npCI(percentageChange(BB.segAmpMid,stm))];
                        burPPC(:,stm,p,cond) = [nanmedian(percentageChange(BB.segPLV,stm)) npCI(percentageChange(BB.segPLV,stm))];
                        burInt(:,stm,p,cond) = [nanmedian(percentageChange(BB.segInterval,stm)) npCI(percentageChange(BB.segInterval,stm))];
                        burRP{1,1,stm,p,cond} = squeeze(BB.segRP{stm}(:,3,4));
                        burRP{1,2,stm,p,cond} = squeeze(BB.segRP{stm}(:,1,4));
                        burRP{2,1,stm,p,cond} = squeeze(BBH.segRP{stm}(:,3,4));
                        burRP{2,2,stm,p,cond} = squeeze(BBH.segRP{stm}(:,1,4));
                        
                        plvStore{1,1,stm,p,cond} = abs(sum(exp(sqrt(-1)*burRP{1,1,stm,p,cond})))./numel(burRP{1,1,stm,p,cond});
                        plvStore{1,2,stm,p,cond} = abs(sum(exp(sqrt(-1)*burRP{1,2,stm,p,cond})))./numel(burRP{1,2,stm,p,cond});
                        plvStore{2,1,stm,p,cond} = abs(sum(exp(sqrt(-1)*burRP{2,1,stm,p,cond})))./numel(burRP{2,1,stm,p,cond});
                        plvStore{2,2,stm,p,cond} = abs(sum(exp(sqrt(-1)*burRP{2,2,stm,p,cond})))./numel(burRP{2,2,stm,p,cond});
                        
                        durStore{1,stm,p} = BB.segDur{stm};
                        ampStore{1,stm,p} = BB.segAmp{stm};
                        ppcStore{1,stm,p} = BB.segPLV{stm};
                        siStore{1,stm,p} = BB.segInterval{stm};
                        
                        durStore{2,stm,p} = BBH.segDur{stm};
                        ampStore{2,stm,p} = BBH.segAmp{stm};
                        ppcStore{2,stm,p} = BBH.segPLV{stm};
                        siStore{2,stm,p} = BBH.segInterval{stm};
                        
                        %                     trajStore{stm,p} = BB.segTraj{stm};
                    else
                        burRate(:,stm,p,cond) = nan;
                        burdur(:,stm,p,cond) = nan;
                        burAmp(:,stm,p,cond) = nan;
                        burAmpMid(:,stm,p,cond) = nan;
                        burPPC(:,stm,p,cond) =nan;
                        burInt(:,stm,p,cond) = nan;
                        
                        durStore{stm,p} = nan;
                        ampStore{stm,p} = nan;
                        ppcStore{stm,p} = nan;
                        siStore{stm,p} =nan;
                        %                     trajStore{stm,p} = BB.segTraj{stm};
                        
                    end
                end
                disp([p,cond,ctype])
            end % Phase shift loop
            
            
            mkdir([Rorg.rootn '\data\CloseLoop_stateDependency'])
            save([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '_' stim_sens '_thresholdFitted.mat'],'powspec_save',...
                'intpow','maxpow','burRate','burdur','burAmp','burAmpMid','burPPC',...
                'durStore','ampStore','ppcStore','siStore',...
                'burInt','phaseShift','NcS','burRP','plvStore'); %,'trajStore')
            
            if cond==1
                save([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '_' stim_sens '_xsims.mat'],'xsimStore','uexsStore')
            end
            
        end % Connection shift loop
    end
    
end

%% Now Plot Results
ctype = 1;
load([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '_' stim_sens '_thresholdFitted.mat'],'powspec_save',...
    'intpow','maxpow','burRate','burdur','burAmp','burAmpMid','burPPC',...
    'durStore','ampStore','ppcStore','siStore',...
    'burInt','phaseShift','NcS','burRP','plvStore'); %,'trajStore')
%% First round of plots
% cmapDisc = brewermap(9,'Set1');
% cmap = brewermap(numel(phaseShift)+4,'Reds');
% cmap = cmap(4:end,:);
cmap = brewermap(12,'Set1');
baseCon = find(NcS==1,1);
figure(1)
phaseShift = rad2deg(phaseShift);
for C =1:3
    if C == 1
        subplot(3,3,1) % Spectra Plot
        titbit = 'M2 Power';
    elseif C == 2
        subplot(3,3,2) % Spectra Plot
        titbit = 'STN Power';
    elseif C == 3
        subplot(3,3,3) % Spectra Plot
        titbit = 'STN/M2 Coherence';
    end
    phsel = 1:2:12;
    ip = 0;
    a = [];
    for i = phsel
        ip = ip+1;
        a(ip) = plot(R.frqz,squeeze(powspec_save(:,C,2,i,baseCon)),'color',cmap(i,:),'LineWidth',2);
        hold on
    end
    plot(R.frqz,squeeze(powspec_save(:,C,1,i,baseCon)),'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    
    xlim([8 34])
    legend(a,sprintfc('%.1f rad.',phaseShift(phsel)))
    xlabel('Frequency (Hz)'); ylabel([titbit])
    title([titbit])
    grid on; axis square
    
    % ARCs
    subplot(3,3,C+3) % Steady State Stats
    % Beta 1
    X = squeeze(intpow(C,1,:,:,baseCon))';
    X = 100.*(X(:,2)-X(:,1))./X(:,1);
    p(1) = plot(phaseShift,X,'color','k','LineWidth',1.5);
    hold on
    s = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');
    % Beta 2
    X = squeeze(intpow(C,2,:,:,baseCon))';
    X = 100.*(X(:,2)-X(:,1))./X(:,1);
    p(2) = plot(phaseShift,X,'color','k','LineWidth',1.5,'LineStyle','--');
    hold on
    s(2) = scatter(phaseShift(phsel),X(phsel),75,cmap(phsel,:),'filled');
    s(2).Marker = 'square';
    
    grid on;axis square
    
    xlabel('Stimulation Phase (radians)'); ylabel('Percentage Change')
    %     legend(p,{'\beta_1 (14-21 Hz) Power','\beta_2 (21-30 Hz) Power'})
    title([titbit ' Response Curve'])
    a = gca;
    a.XTick = rad2deg(([0 pi/2 pi 3*pi/2 2*pi]));
    
    xlim([0 360])
    ylim([-50 150]);
end
set(gcf,'Position',[ 711   418   957   560])

figure(2)
 % High Beta RP
titnames = {'Unstimulated','Max. Suppressing','Max. Amplifying'};
psel = [1 3; 2 7; 2 2];
for i = 1:3
    subplot(2,4,i)
    A = burRP{1,2,psel(i,1),psel(i,2),1};
    p = polarhistogram(A,-pi:pi/8:pi);
    if i == 1
        p.FaceColor = [0 0 0];
    else
        p.FaceColor = cmap(psel(i,2),:);
    end
    plv(i) = abs(sum(exp(sqrt(-1)*A)))./numel(A);
    title([titnames{i} sprintf('\n PLV = %0.3f',plv(i))] )
end
subplot(4,4,4)
bar(plv)

 % High Beta RP
psel = [1 3; 2 2; 2 7];
for i = 1:3
    subplot(2,4,i+4)
    A = burRP{2,2,psel(i,1),psel(i,2),1};
    p = polarhistogram(A,-pi:pi/8:pi);
    if i == 1
        p.FaceColor = [0 0 0];
    else
        p.FaceColor = cmap(psel(i,2),:);
    end
    plv(i) = abs(sum(exp(sqrt(-1)*A)))./numel(A);
    title([titnames{i} sprintf('\n PLV = %0.3f',plv(i))] )
end
subplot(2,4,8)
bar(plv)

figure(1)
ip = 0;
for  i = [2 3 4]
    ip = ip +1;
    if i == 1
        X = squeeze(burRate(1,:,:,baseCon))';
        Y = nan(size(X));
        titname = 'Burst Rate';
        ylab = 'Burst Probability (sec-1)';
        ls = '-';
        rlz = [-50 50]; %[1.25 1.75].*10^-6;
    elseif i == 2
        X = squeeze(burInt(1,:,:,baseCon))';
        Y = squeeze(burInt(2,:,:,baseCon))';
        titname = 'Inter-Burst Interval';
        ylab = 'IBI (ms)';
        rlz = [-75 75]; %[1.25 1.75].*10^-6;
    elseif i == 3
        X = squeeze(burAmp(1,:,:,baseCon))';
        Y = squeeze(burAmp(2,:,:,baseCon))';
        titname = 'Burst Amplitude';
        ylab = 'Peak Amplitude (a.u.)';
        rlz = [-15 15]; %[1.25 1.75].*10^-6;
    elseif i == 4
        X = squeeze(burdur(1,:,:,baseCon))';
        Y = squeeze(burdur(2,:,:,baseCon))';
        titname = 'Burst Duration';
        ylab = 'Duration (ms)';
        rlz = [-75 125]; %[1.25 1.75].*10^-6;
    elseif i == 5
        X = squeeze(burPPC(1,:,:,baseCon))';
        Y = squeeze(burPPC(2,:,:,baseCon))';
        titname = 'Burst STN/M2 Synchronization';
        ylab = 'PPC Magnitude';
        rlz = [-50 50]; %[1.25 1.75].*10^-6;
    end
    %
    Z = X(:,1);
    X = X(:,2);
    
    subplot(3,3,6+ip)
    [l b] = boundedline((phaseShift),X,Y(:,2)./2);
    l.Color = 'k';
    l.LineWidth = 1.5;
    b.FaceColor = 'k';
    b.FaceAlpha = 0.2;
    hold on
    s = scatter((phaseShift(phsel)),X(phsel),75,cmap(phsel,:),'filled');
    
    plot(phaseShift,Z,'LineWidth',1,'Color','k','LineStyle','-'); % Unstimulated median
    plot(phaseShift,Z-Y(:,1)./2,'LineWidth',1,'Color','k','LineStyle','--'); % SEM
    plot(phaseShift,Z+Y(:,1)./2,'LineWidth',1,'Color','k','LineStyle','--')
    
    xlabel('Stimulation Phase')
    ylabel(ylab)
    title(titname)
    a = gca;
    a.XTick = rad2deg([0 pi/2 pi 3*pi/2 2*pi]);
    
    xlim([0 360])
    grid on; axis square
    ylim(rlz);
    % rlim(rlz)
end
set(gcf,'Position',[451   13  908  764])

%% Statistics for Bursts
X = squeeze(intpow(2,1,:,:,1))';

[a pSup] = min(X(:,2));
[a pAmp] = max(X(:,2));

[h p ci stats] = ttest2(durStore{1,pSup},durStore{2,pSup});
durSV = [nanmean(durStore{1,pSup})-nanmean(durStore{2,pSup}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(ampStore{1,pSup},ampStore{2,pSup});
ampSV = [nanmean(ampStore{1,pSup})-nanmean(ampStore{2,pSup}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(log(siStore{1,pSup}(5:end)),log(siStore{2,pSup}(5:end)));
siSV = [nanmean(siStore{1,pSup})-nanmean(siStore{2,pAmp}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(log(ppcStore{1,pSup}(3:end)),log(ppcStore{2,pSup}(3:end)));
ppcSV = [nanmean(ppcStore{1,pSup}(3:end))-nanmean(ppcStore{2,pSup}(3:end)) diff(ci)./2 stats.df stats.tstat p 0.05/5]


%% Trajectory Plots
% figure
% subplot(1,3,1)
% plot(nanmean(squeeze(trajStore{2,pAmp}(4,:,:)),2)); hold on; plot(nanmean(squeeze(trajStore{2,pSup}(4,:,:)),2))
% legend({'Amplifying','Supressing'});
% subplot(1,3,2)
% plot(nanmean(squeeze(trajStore{1,pAmp}(4,:,:)),2)); hold on; plot(nanmean(squeeze(trajStore{2,pAmp}(4,:,:)),2))
% legend({'Base','Amplifying'});
% subplot(1,3,3)
% plot(nanmean(squeeze(trajStore{1,pSup}(4,:,:)),2)); hold on; plot(nanmean(squeeze(trajStore{2,pSup}(4,:,:)),2))
% legend({'Base','Supressing'});
%
function ci = npCI(X)
N = numel(X);
lb = (N/2) - ((1.96*sqrt(N))/2);
ub = 1+ (N/2) + ((1.96*sqrt(N))/2);

X = sort(X);
ub = fix(ub); lb = fix(lb);
if (ub>0) & (lb>0)
    ub = X(fix(ub));
    lb = X(fix(lb));
    ci = ub-lb;
else
    ci = nan;
end

function prcD = percentageChange(PLR,stm)
prcD = ((PLR{stm}-nanmedian(PLR{1}))./nanmedian(PLR{1}))*100; % Make Amp % Change

% % % Burst Duration
% % X = burdur';
% % X = 100.*(X(:,2)-X(:,1))./X(:,1);
% % p(1) = plot(phaseShift,X,'color',cmap(6,:),'LineWidth',2);
% % hold on
% % s(1) = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');
% %
% % % Burst Amplitude
% % X = burAmp';
% % X = 100.*(X(:,2)-X(:,1))./X(:,1);
% % p(2) = plot(phaseShift,X,'color',cmap(12,:),'LineWidth',2);
% % hold on
% % s(2) = scatter(phaseShift(phsel),X(phsel),65,cmap(phsel,:),'filled');
% % s(2).Marker = 'diamond';
% %
% % % Burst Phase Locking
% % X = burPPC';
% % X = 100.*(X(:,2)-X(:,1))./X(:,1);
% % p(3) = plot(phaseShift,X,'color',cmap(18,:),'LineWidth',2);
% % hold on
% % s(3) = scatter(phaseShift(phsel),X(phsel),65,cmap(phsel,:),'filled');
% % s(3).Marker = 'square';
% %

%% For repetitive pulses
% pulseAmp = 0.05;
% pulseCycle = (1/50).*fsamp;
% pulseDuty = (0.005*fsamp);
% % If rhythmic Stim
% %         pulseKern = zeros(pulseWid,1);
% %         ps = 1;
% %         while ps <= pulseWid
% %         pulseKern(ps:ps+pulseDuty) = pulseAmp;
% %         ps = ps + floor(pulseCycle+ 5.*randn);
% %         end
% pulseKern = ones(pulseWid,1);
