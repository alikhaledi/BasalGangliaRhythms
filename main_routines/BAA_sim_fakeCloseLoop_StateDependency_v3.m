function [R] = BAA_sim_fakeCloseLoop_StateDependency_v3(Rorg,modID,simtime,fresh)
% Load in model data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v3(R,modID,simtime);]
% or load from preload
load([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
R.filepathn =  Rorg.filepathn;
R.rootn =  Rorg.rootn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')
% pause(1)
for SScomb = 1:2
    %% Define stimulation conditions
    if SScomb == 1
        % Stimualating M2
        senssite = 4; % STN
        stimsite = 1; % M2
        stim_sens = 'stimM2_sensSTN';
    elseif SScomb == 2
        % Stimulating  STN
        senssite = 1; % M2
        stimsite = 4; % STN
        stim_sens = 'stimSTN_sensM2';
    end
    
    R.IntP.phaseStim.sensStm = [senssite stimsite];
    
    P = permMod{1}.par_rep{1};
    % Simulation Coniditions
    R.obs.csd.df = 0.5;
    R = setSimTime(R,simtime);
    R.obs.trans.norm = 0; % No normalization of spectra
    R.obs.brn = R.obs.brn;
    R.IntP.intFx = @spm_fx_compile_120319_stim;
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    tplot = 0;
    % Give all timeseries the same input - makes comparable
    rng(5453)
    uc = innovate_timeseries(R,m);
    uc{1} = uc{1}.*sqrt(R.IntP.dt);
    
    if fresh
        for ctype = 1:2
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
            
            for cond = 1;%:numel(NcS)
                
                % Get Base Parameters
                Pbase = P;
                
                % Modulate connectivity (state dependency)
                if ctype ==1
                    Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*NcS(cond)); %
                elseif ctype == 2
                    Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*NcS(cond)); %
                end
                
                %% Setup stim parameterss
                switch stim_sens
                    case 'stimM2_sensSTN'
                        R = typeIstimPars_v3(R);
                    case 'stimSTN_sensM2'
                        %                     R = typeIIstimPars(R);
                        R = typeIstimPars_v3(R);
                end
                
                %% Simulate Base Model
                uc_ip{1} = uc;
                R.IntP.phaseStim.switch = 0 ;
                [~,~,feat_sim{1},xsim_gl,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
                
                %% Work out the threshold
                [~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);
                eps(cond) =  R.IntP.phaseStim.eps;
                
                for p = 1:numel(phaseShift) %[1 10] %
                    %% Resimulate with Phase-Locked Input
                    R.IntP.phaseStim.switch = 1;
                    
                    % Modulate the phase
                    R.IntP.phaseStim.phaseshift = phaseShift(p);
                    
                    % Simulate with Stimulation
                    [~,~,feat_sim{2},~,xsim_ip{2},~,Rout]  = computeSimData(R,m,uc_ip{1},Pbase,0);
                    load([R.rootn 'data\rat_InDirect_ModelComp\phaseStimSave\stim_tmp'],'uexs')
                    
                    pU = uexs(R.IntP.phaseStim.sensStm(2),round(R.obs.brn*(1/R.IntP.dt))+1:end);
                    indsS = SplitVec(find(abs(pU)>0),'consecutive'); % Split up data based upon the target threshold
                    indsN = SplitVec(find(abs(pU)==0),'consecutive'); % Split up data based upon the target threshold
                    
                    %                     for c = 1:numel(indsS)-1
                    %                     indsS{c} = [indsS{c} indsS{c}(end):indsS{c}(end)+300];
                    %                     end
                    %                     indsS{end} = [];
                    
                    XS{1} = xsim_ip{1}{1}([1 4],[indsS{:}]); % No stim
                    XS{2} = xsim_ip{2}{1}([1 4],[indsS{:}]); % with stim
                    
                    tplot = 0;
                    if tplot == 1
                        pU = uexs(1,round(R.obs.brn*(1/R.IntP.dt))+1:end);
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
                        xlim([30.9 31.9])
                        R.plot.outFeatFx({feat_sim{1}},{feat_sim{2}},R.data.feat_xscale,R,1,[])
                        
                    end
                    
                    for stm = 1:2
                        [F,Hz] = pwelch(XS{stm}',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
                        [C,Hz] = mscohere(XS{stm}(1,:)',XS{stm}(2,:)',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
                        
                        spec = [F C];
                        
                        powspec_save(:,:,stm,p,cond) = spec;
                        intpow(:,1,stm,p,cond) = sum(spec(Hz>14 & Hz<=21,:))./numel(find(Hz>14 & Hz<=30));
                        maxpow(:,1,stm,p,cond) = max(spec(Hz>14 & Hz<=21,:));
                        intpow(:,2,stm,p,cond) = sum(spec(Hz>21 & Hz<=30,:))./numel(find(Hz>14 & Hz<=30));
                        maxpow(:,2,stm,p,cond) = max(spec(Hz>21 & Hz<=30,:));
                        intpow(:,3,stm,p,cond) = sum(spec(Hz>14 & Hz<=30,:))./numel(find(Hz>14 & Hz<=30));
                        maxpow(:,3,stm,p,cond) = max(spec(Hz>14 & Hz<=30,:));
                        
                        xsimStore{stm,p} = xsim_ip{stm}{1};
                        uexsStore{stm,p} = uexs;
                    end
                    
                    
                    disp([p,cond,ctype])
                end % Phase shift loop
                
                
                mkdir([Rorg.rootn '\data\CloseLoop_stateDependency'])
                save([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_v3_save_' num2str(ctype) '_' stim_sens '_thresholdFitted.mat'],'powspec_save',...
                    'intpow','maxpow','burRate','burdur','burAmp','burAmpMid','burPPC',...
                    'durStore','ampStore','ppcStore','siStore',...
                    'burInt','phaseShift','NcS','Hz'); %,'burRP','plvStore'); %,'trajStore')
                
                if cond==1
                    save([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_v3_save_' num2str(ctype) '_' stim_sens '_xsims.mat'],'xsimStore','uexsStore')
                end
                
            end
        end
    end
    ctype = 1;
    load([Rorg.rootn '\data\CloseLoop_stateDependency\CloseLoop_stateDependency_v3_save_' num2str(ctype) '_' stim_sens '_thresholdFitted.mat'],'powspec_save',...
        'intpow','maxpow','burRate','burdur','burAmp','burAmpMid','burPPC',...
        'durStore','ampStore','ppcStore','siStore',...
        'burInt','phaseShift','NcS','Hz'); %,'burRP','plvStore'); %,'trajStore')
    
    figure(SScomb)
    closedLoopSpectralPlot(R,phaseShift,NcS,intpow,powspec_save,Hz)
    
end
% for ch = 1:2
%     for stm = 1:2
%         subplot(2,2,ch)
%         plot(squeeze(peakLBPow(ch,stm,:,1)))
%         hold on
%
%         subplot(2,2,ch+2)
%         plot(squeeze(peakHBPow(ch,stm,:,1)))
%         hold on
%
%     end
% end
