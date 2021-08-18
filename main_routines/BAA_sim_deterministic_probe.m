function [R] = BAA_sim_deterministic_probe(Rorg,modID,simtime,fresh)
% Load in model data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v3(R,modID,simtime);]

% Stimualating M2
stimsite = 1; % M2
senssite = 4; % STN
stim_sens = 'stimM2_sensSTN';
%
% Stimulating  STN
% stimsite = 4; % STN
% senssite = 1; % M2
% stim_sens = 'stimSTN_sensM2'
%


% or load from preload
load([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
warning('Loading Preloaded model, cant change simtime or model choice!!!')
pause(1)

P = permMod{1}.par_rep{1};
% Simulation Coniditions
R.obs.csd.df = 0.5;
R = setSimTime(R,simtime);
R.obs.trans.norm = 1; % No normalization of spectra
R.IntP.intFx = @spm_fx_compile_120319;
phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12
tplot = 0;
% Give all timeseries the same input - makes comparable
R.IntP.Utype = 'zero';
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
fsamp = 1/R.IntP.dt;
% uc{1}(:,1) = 0.125;
uc{1}((R.obs.brn/R.IntP.dt):(2/R.IntP.dt):end,1) = 0.025;

if fresh
    for ctype = 1:2
        uc_ip{1} = uc;
        Pbase = P;
        %% Simulate Base Model
        [~,~,feat_sim{1},xsim_ip_gl{1},xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0,1);
        outstates = [1:2:7 9 11 13 15 17];
        
        init = find(abs(sum(xsim_ip_gl{1}{1},1))>0,1,'first');
        window = init:init+(0.5/R.IntP.dt);
        
        figure(100)
        % plot(xsim_ip_gl{1}{1}([11],window)',xsim_ip_gl{1}{1}([13],window)')
        plot(xsim_ip_gl{1}{1}([11],window)',xsim_ip_gl{1}{1}([13],window)');
        
        %         Pbase.int{1}.G(6) = -32;
        %         Pbase.int{1}.G(5) = -32;
        for p = 1:numel(phaseShift) %[1 10] %
            % Now find bursts that will be used to parameterize
            % stimulation
            R.condname = {'1'};
            [R,BB] = compute_BetaBursts_Simulated(R,xsim_ip{1},0);
            if cond == 1 && ctype == 1
                R.BB.thresh_prctile = 75;% o85; tl 80
                R.BB.threshold_type = 'baseModelThresh';
                RBB = compute_BurstThreshold(R,BB,1,0);
            end
            
            BB.epsAmpfull = RBB.epsAmpfull;
            BB.epsAmp = RBB.epsAmp;
            BB.epsPLV = RBB.epsPLV;
            
            R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
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
                    pulse_Phi = BB.Phi{1}(senssite,pulseInds);
                    pulseKern = sin(pulse_Phi+(phaseShift(p))); %
                    pU(pulseInds) = pulseKern;
                end
            end
            pU = (1/4.*std(uc{1}(:,stimsite))).*pU; %.*pulseAmp;
            uc_ip{2} =  uc_ip{1};
            uc_ip{2}{100} = zeros(size(uc_ip{1}{1}));
            uc_ip{2}{100}(:,stimsite) =  pU'; % Give it a cortical pulse
            
            % Simulate with Stimulation
            [~,~,feat_sim{2},~,xsim_ip{2},~,Rout]  = computeSimData(R,m,uc_ip{2},Pbase,0);
            
            if tplot == 1
                % Optional Plots for TimeSeries
                figure(100)
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
            R.condname = {'1','2'};
            [R,BB] = compute_BetaBursts_Simulated(R,{xsim_ip{1}{1} xsim_ip{2}{1}});
            %                 R.BB.thresh_prctile = 75;% o85; tl 80
            %                 R.BB.threshold_type = 'baseModelThresh';
            %                 BB = compute_BurstThreshold(R,BB,1,0);
            
            BB.epsAmpfull = RBB.epsAmpfull;
            BB.epsAmp = RBB.epsAmp;
            BB.epsPLV = RBB.epsPLV;
            
            R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
            BB.plot.durlogflag = 0;
            R.BB.pairInd = [1 4]; % KEEP LOOKING AT STN BURSTS!
            BB = defineBetaEvents(R,BB);
            
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
                
                if numel(BB.segAmp)>2
                    burRate(:,stm,p,cond) = percentageChange(BB.segRate,stm);
                    burdur(:,stm,p,cond) = [nanmedian(percentageChange(BB.segDur,stm)) npCI(percentageChange(BB.segDur,stm))];
                    burAmp(:,stm,p,cond) = [nanmedian(percentageChange(BB.segAmp,stm)) npCI(percentageChange(BB.segAmp,stm))];
                    burAmpMid(:,stm,p,cond) = [nanmedian(percentageChange(BB.segAmpMid,stm)) npCI(percentageChange(BB.segAmpMid,stm))];
                    burPPC(:,stm,p,cond) = [nanmedian(percentageChange(BB.segPLV,stm)) npCI(percentageChange(BB.segPLV,stm))];
                    burInt(:,stm,p,cond) = [nanmedian(percentageChange(BB.segInterval,stm)) npCI(percentageChange(BB.segInterval,stm))];
                    
                    durStore{stm,p} = BB.segDur{stm};
                    ampStore{stm,p} = BB.segAmp{stm};
                    ppcStore{stm,p} = BB.segPLV{stm};
                    siStore{stm,p} = BB.segInterval{stm};
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
    end % Connection shift loop
    mkdir([Rorg.rootn '\data\CloseLoop_stateDependency'])
    save([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_save_' num2str(ctype) '_' stim_sens '_thresholdFitted.mat'],'powspec_save',...
        'intpow','maxpow','burRate','burdur','burAmp','burAmpMid','burPPC',...
        'durStore','ampStore','ppcStore','siStore',...
        'burInt','phaseShift','NcS'); %,'trajStore')
    
end



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
        

