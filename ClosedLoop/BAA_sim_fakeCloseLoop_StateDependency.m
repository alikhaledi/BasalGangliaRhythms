function [R] = BAA_sim_fakeCloseLoop_StateDependency(R,modID,simtime,fresh)
% Load in model data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v3(R,modID,simtime);
stimsite = 4; % STN
senssite = 1; % M2

% or load from preload
load([R.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
warning('Loading Preloaded model, cant change simtime or model choice!!!')
pause(1)

P = permMod{1}.par_rep{1};
% Simulation Coniditions
R.obs.csd.df = 0.5;
R = setSimTime(R,simtime);
R.obs.trans.norm = 0; % No normalization of spectra
NcS = 18;
phaseShift = linspace(0,2.*pi,14); % List of phases to be tested

tplot = 0;
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
fsamp = 1/R.IntP.dt;

if fresh
    for ctype =1:2
        % temp!
%         load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '.mat'],'powspec_save','intpow','maxpow','burRate','burdur','burAmp','burPPC',...
%             'durStore','ampStore','ppcStore','siStore',...
%             'burInt','phaseShift','conStren')
        
        if ctype ==1
            conStren = [1 linspace(0.05,1.5,NcS)]; % [0.1 1 1.15]; STN-> GPe
        elseif ctype == 2
            conStren = [1 linspace(0.05,2.5,NcS)]; % [0.1 1 1.15]; M2 -> STN
        end
        for cond = 1;%:NcS+1
            
            % Get Base Parameters
            Pbase = P;
            R.obs.brn = 0; % temporarily!
            
            % Modulate connectivity (state dependency)
            if ctype ==1
                Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*conStren(cond)); %
            elseif ctype == 2
                Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*conStren(cond)); %
            end
            %     intpow = []; maxpow = [];
            for p = 1:numel(phaseShift) %[1 10] %
                
                % Initialize variables
                uc_ip = {}; feat_sim = {}; xsim_ip = {};
                uc_ip{1} = uc;
                
                %% Simulate Base Model
                [~,~,feat_sim{1},~,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
                
                % Now find bursts that will be used to parameterize
                % stimulation
                R.condname = {'1'};
                [R,BB] = compute_BetaBursts_Simulated(R,xsim_ip{1},0);
                R.BB.thresh_prctile = 75;% o85; tl 80
                BB = compute_BurstThreshold(R,BB,1,0);
                R.BB.minBBlength = 0.5; %o1 tl 1.5; %  Minimum burst period- cycles
                BB.plot.durlogflag = 0;
                R.BB.pairInd = [1 senssite]; % Use cortical burst sensing
                BB = defineBetaEvents(R,BB);
                
                
                %% Resimulate with Phase-Locked Input
                % Now decide the timing of the intervention!
                pulseWid = (0.3.*fsamp);
                pulseDelay = (0.*fsamp);
                
                % Setup the Stimulation Input
                pU = zeros(size(R.IntP.tvec));
                pulseStart = [];
                for seg = 1:numel(BB.segInds{1})
                    pulseStart(seg) = BB.segInds{1}{seg}(1) + pulseDelay;
                    pulseInds = pulseStart(seg):pulseStart(seg)+pulseWid-1;
                    if pulseInds(end)<=size(BB.Phi{1},2) % Ensure always within sample range
                        pulse_Phi = BB.Phi{1}(4,pulseInds);
                        pulseKern = sin(pulse_Phi+(phaseShift(p))); %
                        pU(pulseInds) = pulseKern;
                    end
                end
                pU = (0.25.*std(uc{1}(:,stimsite))).*pU; %.*pulseAmp;
                uc_ip{2} =  uc_ip{1};
                uc_ip{2}{1}(:,stimsite) = uc_ip{2}{1}(:,stimsite) + pU'; % Give it a cortical pulse
                
                % Simulate with Stimulation
                [~,~,feat_sim{2},~,xsim_ip{2},~,Rout]  = computeSimData(R,m,uc_ip{2},Pbase,0);
                
                if tplot == 1
                    % Optional Plots for TimeSeries
                    figure(100)
                    a(1) = subplot(2,1,1);
                    plot(Rout.IntP.tvec_obs,xsim_ip{1}{1}(1,2:end));
                    hold on
                    plot(Rout.IntP.tvec_obs,xsim_ip{2}{1}(1,2:end));
                    
                    a(2) = subplot(2,1,2);
                    plot(Rout.IntP.tvec_obs,xsim_ip{1}{1}(4,2:end)+5e-6);
                    hold on
                    plot(Rout.IntP.tvec_obs,xsim_ip{2}{1}(4,2:end)+5e-6);
                    linkaxes(a,'x')
                    xlim([7 8])
                    R.plot.outFeatFx({feat_sim{1}},{feat_sim{2}},R.data.feat_xscale,R,1,[])
                end
                
                % Re-compute bursts (simulated data)
                R.condname = {'1','2'};
                [R,BB] = compute_BetaBursts_Simulated(R,{xsim_ip{1}{1} xsim_ip{2}{1}});
                R.BB.thresh_prctile = 75;% o85; tl 80
                BB = compute_BurstThreshold(R,BB,1,0);
                R.BB.minBBlength = 0.5; %o1 tl 1.5; %  Minimum burst period- cycles
                BB.plot.durlogflag = 0;
                R.BB.pairInd = [1 senssite]; % Use cortical burst sensing
                BB = defineBetaEvents(R,BB);
                
                %% Now get statistics of the simulations (spectra/burst features)
                
                for stm = 1:2
                    spec = [squeeze(feat_sim{stm}(1,1,1,1,:)) squeeze(feat_sim{stm}(1,4,4,1,:)) squeeze(feat_sim{stm}(1,4,1,4,:))];
                    powspec_save(:,:,stm,p,cond) = spec;
                    intpow(:,1,stm,p,cond) = sum(spec(R.frqz>14 & R.frqz<=21,:));
                    maxpow(:,1,stm,p,cond) = max(spec(R.frqz>14 & R.frqz<=21,:));
                    intpow(:,2,stm,p,cond) = sum(spec(R.frqz>21 & R.frqz<=30,:));
                    maxpow(:,2,stm,p,cond) = max(spec(R.frqz>21 & R.frqz<=30,:));
                    intpow(:,3,stm,p,cond) = sum(spec(R.frqz>21 & R.frqz<=30,:));
                    maxpow(:,3,stm,p,cond) = max(spec(R.frqz>21 & R.frqz<=30,:));
                    
                    burRate(:,stm,p,cond) = percentageChange(BB.segRate,stm);
                    burdur(:,stm,p,cond) = [nanmedian(percentageChange(BB.segDur,stm)) npCI(percentageChange(BB.segDur,stm))];
                    burAmp(:,stm,p,cond) = [nanmedian(percentageChange(BB.segAmp,stm)) npCI(percentageChange(BB.segAmp,stm))];
                    burPPC(:,stm,p,cond) = [nanmedian(percentageChange(BB.segPLV,stm)) npCI(percentageChange(BB.segPLV,stm))];
                    burInt(:,stm,p,cond) = [nanmedian(percentageChange(BB.segInterval,stm)) npCI(percentageChange(BB.segInterval,stm))];
                    
                    durStore{stm,p} = BB.segDur{stm};
                    ampStore{stm,p} = BB.segAmp{stm};
                    ppcStore{stm,p} = BB.segPLV{stm};
                    siStore{stm,p} = BB.segInterval{stm};
%                     trajStore{stm,p} = BB.segTraj{stm};
                end
                
                
                
                disp([p,cond,ctype])
            end % Phase shift loop
        end % Connection shift loop
        mkdir([R.rootn '\data\CloseLoop_stateDependency'])
        save([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '.mat'],'powspec_save','intpow','maxpow','burRate','burdur','burAmp','burPPC',...
            'durStore','ampStore','ppcStore','siStore',...
            'burInt','phaseShift','conStren','trajStore')
    end
    
end

%% Now Plot Results
ctype = 1;
load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '.mat'],'powspec_save','intpow','maxpow','burRate','burdur','burAmp','burPPC',...
    'durStore','ampStore','ppcStore','siStore',...
    'burInt','phaseShift','conStren','trajStore')
%% First round of plots
baseCon = 1;
cmapDisc = brewermap(9,'Set1');
cmap = brewermap(numel(phaseShift)+4,'Reds');
cmap = cmap(4:end,:);
% phaseShift = rad2deg(phaseShift);
for C =1:3
    if C == 1
        subplot(2,3,1) % Spectra Plot
        titbit = 'M2 Power';
    elseif C == 2
        subplot(2,3,2) % Spectra Plot
        titbit = 'STN Power';
    elseif C == 3
        subplot(2,3,3) % Spectra Plot
        titbit = 'STN/M2 Coherence';
    end
    phsel = 1:3:14;
    ip = 0;
    a = [];
    for i = phsel
        ip = ip+1;
        a(ip) = plot(R.frqz,squeeze(powspec_save(:,C,2,i,baseCon)),'color',cmap(i,:),'LineWidth',2);
        hold on
    end
    plot(R.frqz,squeeze(powspec_save(:,C,1,i,baseCon)),'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    xlim([8 34])
    %     legend(a,sprintfc('%.1f rad.',phaseShift(phsel)))
    xlabel('Frequency (Hz)'); ylabel([titbit])
    title([titbit])
    grid on
    
    % ARCs
    subplot(2,3,C+3) % Steady State Stats
    % Beta 1
    X = squeeze(intpow(C,1,:,:,1))';
    X = 100.*(X(:,2)-X(:,1))./X(:,1);
    p(1) = plot(phaseShift,X,'color','k','LineWidth',1.5);
    hold on
    s = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');
    % Beta 2
    X = squeeze(intpow(C,2,:,:,1))';
    X = 100.*(X(:,2)-X(:,1))./X(:,1);
    p(2) = plot(phaseShift,X,'color','k','LineWidth',1.5,'LineStyle','--');
    hold on
    s(2) = scatter(phaseShift(phsel),X(phsel),75,cmap(phsel,:),'filled');
    s(2).Marker = 'square';
    
    grid on
    
    xlabel('Stimulation Phase (radians)'); ylabel('Percentage Change')
    %     legend(p,{'\beta_1 (14-21 Hz) Power','\beta_2 (21-30 Hz) Power'})
    title([titbit ' Response Curve'])
    a = gca;
    a.XTick = ([0 pi/2 pi 3*pi/2 2*pi]);
    
    xlim([0 2*pi])
    ylim([-50 150]);
end
set(gcf,'Position',[ 711   418   957   560])

figure
for  i = 1:5
    if i == 1
        X = squeeze(burRate(1,:,:,1))';
        Y = nan(size(X));
        titname = 'Burst Rate';
        ylab = 'Burst Probability (sec-1)';
        ls = '-';
        rlz = [-50 50]; %[1.25 1.75].*10^-6;
    elseif i == 2
        X = squeeze(burInt(1,:,:,1))';
        Y = squeeze(burInt(2,:,:,1))';
        titname = 'Inter-Burst Interval';
        ylab = 'IBI (ms)';
        rlz = [-50 50]; %[1.25 1.75].*10^-6;
    elseif i == 3
        X = squeeze(burAmp(1,:,:,1))';
        Y = squeeze(burAmp(2,:,:,1))';
        titname = 'Burst Amplitude';
        ylab = 'Peak Amplitude (a.u.)';
        rlz = [-15 15]; %[1.25 1.75].*10^-6;
    elseif i == 4
        X = squeeze(burdur(1,:,:,1))';
        Y = squeeze(burdur(2,:,:,1))';
        titname = 'Burst Duration';
        ylab = 'Duration (ms)';
        rlz = [-50 50]; %[1.25 1.75].*10^-6;
    elseif i == 5
        X = squeeze(burPPC(1,:,:,1))';
        Y = squeeze(burPPC(2,:,:,1))';
        titname = 'Burst STN/M2 Synchronization';
        ylab = 'PPC Magnitude';
        rlz = [-50 50]; %[1.25 1.75].*10^-6;
    end
    %
    Z = X(:,1);
    X = X(:,2);
    
    subplot(2,3,i+1)
    [l b] = boundedline((phaseShift),X,Y(:,2)./2);
    b.FaceAlpha = 0.8;
    hold on
    s = scatter((phaseShift(phsel)),X(phsel),75,cmap(phsel,:),'filled');
    
    plot(phaseShift,Z,'LineWidth',1,'Color','k','LineStyle','-'); % Unstimulated median
    plot(phaseShift,Z-Y(:,1)./2,'LineWidth',1,'Color','k','LineStyle','--'); % SEM
    plot(phaseShift,Z+Y(:,1)./2,'LineWidth',1,'Color','k','LineStyle','--')
    
    xlabel('Stimulation Phase')
    ylabel(ylab)
    title(titname)
    a = gca;
    a.XTick = ([0 pi/2 pi 3*pi/2 2*pi]);
    
    xlim([0 2*pi])
    grid on
    ylim(rlz);
    % rlim(rlz)
end
set(gcf,'Position',[557         238        1235         740])

%% Statistics for Bursts
load([R.rootn '\data\CloseLoop_stateDependency\CloseLoop_burstStats_save_' num2str(ctype) '.mat'],'durStore','ampStore','ppcStore','siStore')
[a pSup] = min(squeeze(burAmp(1,2,:,1)));
[a pAmp] = max(squeeze(burAmp(1,2,:,1)));

[h p ci stats] = ttest2(durStore{1,pSup},durStore{2,pSup});
durSV = [nanmean(durStore{1,pSup})-nanmean(durStore{2,pSup}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(ampStore{1,pSup},ampStore{2,pSup});
ampSV = [nanmean(ampStore{1,pSup})-nanmean(ampStore{2,pSup}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(log(siStore{1,pSup}(5:end)),log(siStore{2,pSup}(5:end)));
siSV = [nanmean(siStore{1,pAmp})-nanmean(siStore{2,pAmp}) diff(ci)./2 stats.df stats.tstat p 0.05/5]

[h p ci stats] = ttest2(log(ppcStore{1,pSup}(3:end)),log(ppcStore{2,pSup}(3:end)));
ppcSV = [nanmean(ppcStore{1,pAmp})-nanmean(ppcStore{2,pAmp}) diff(ci)./2 stats.df stats.tstat p 0.05/5]


%% Trajectory Plots
figure
subplot(1,3,1)
plot(nanmean(squeeze(trajStore{2,pAmp}(4,:,:)),2)); hold on; plot(nanmean(squeeze(trajStore{2,pSup}(4,:,:)),2))
legend({'Amplifying','Supressing'});
subplot(1,3,2)
plot(nanmean(squeeze(trajStore{1,pAmp}(4,:,:)),2)); hold on; plot(nanmean(squeeze(trajStore{2,pAmp}(4,:,:)),2))
legend({'Base','Amplifying'});
subplot(1,3,3)
plot(nanmean(squeeze(trajStore{1,pSup}(4,:,:)),2)); hold on; plot(nanmean(squeeze(trajStore{2,pSup}(4,:,:)),2))
legend({'Base','Supressing'});

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
