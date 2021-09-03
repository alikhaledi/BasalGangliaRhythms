function [R] = BAA_sim_phaseLockedStim_SICompEps(Rorg)
close all
ip = 0;
% Comopute simulations by sweeping across data
% [R,m,permMod] = getSimModelData_v3(Rorg,modID,simtime);
% mkdir([Rorg.rootn 'data\ModelFit\'])
% save([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')

% OR Load it in:
load([Rorg.rootn 'data\modelfit\SimModelData_M10.mat'],'R','m','permMod')
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')

EpsVec = [0:15:100];

for EpsN = 1:numel(EpsVec)
    SScomb =1
    %% Define stimulation conditions
    % Stimualating M2
    senssite = 4; % STN
    stimsite = 1; % M2
    stimFx = @zeroCrossingPhaseStim_v3;
    stim_sens = 'stimM2_sensSTN';
    phflag = 1;
    % Setup stim parameterss
    R = typeIstimPars_v3(R);
    R.IntP.phaseStim.sensStm = [senssite stimsite];
    R.IntP.phaseStim.epsthresh = EpsVec(EpsN); % Overwrite
    R.IntP.phaseStim.stimGap = 0.2;
    R.IntP.phaseStim.stimperiod =1;
    R.IntP.phaseStim.stimAmp = 1/2;
    %% Connection Sets
    ck_1(1,:) = [1 logspace(-1,log10(5),34)];
    ck_1(2,:) = [1 logspace(-1,log10(1.90),34)];

    %     ck_1(1,:) = [1 ck_1_org(1,[2 31])];
    %     ck_1(3,:) = [1 ck_1_org(3,[2 31])];

    %     % Simulation Coniditions
    R.obs.csd.df = 0.5;
    R = setSimTime(R,32);

    % Trans Options
    %     R.obs.SimOrd = 10;
    R.obs.trans.norm = 0;
    R.obs.gainmeth = {};

    % Observe Middle layers
    R.obs.outstates(1) = 3; % change to middle layer
    m.outstates{1} = [0 0 1 0 0 0 0 0];

    % Give all timeseries the same input - makes comparable
    rng(5453)
    m.uset.p.scale = m.uset.p.scale;
    uc = innovate_timeseries(R,m);
    uc{1} = uc{1}.*sqrt(R.IntP.dt);
    XBase = permMod{1}.par_rep{1};

    % Phase To be Tested
    R.IntP.intFx = @spm_fx_compile_120319_stim;

    if phflag
        phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
        phaseShift = phaseShift(1:12); %12
    else
        phaseShift = 0;
    end

    %% Loop through Connections
    CON = 1;
    state = 1
    %% Setup Base Model
    Pbase = XBase;

    % Simulate Base Model
    uc_ip{1} = uc;
    R.IntP.phaseStim.switch = 0 ;
    R.IntP.phaseStim.phaseshift = 0;
    R.frqz = 2:0.2:150;
    R.IntP.compFx = @nullComp;
    [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1},~,Rout] = computeSimData(R,m,uc_ip{1},Pbase,0);
    % Work out the threshold
    R.IntP.phaseStim.eps = 0;
    [~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);

    %% Do Stimulation
    R.IntP.phaseStim.switch = 1;
    m = m; % initialise for parfor
    xsim_ip_stim = cell(1,12); feat_sim_stim = cell(1,12); pU = cell(1,12);
    %             parfor p = 1:numel(phaseShift)
    parfor p = 1:numel(phaseShift)
        Rpar = R;
        % Modulate the phase
        Rpar.IntP.phaseStim.phaseshift = phaseShift(p);
        Rpar.IntP.phaseStim.stimFx = stimFx;
        % Simulate with Stimulation
        [~,~,feat_sim_stim{p},~,xsim_ip_stim{p},~,Rpar]  = computeSimData(Rpar,m,uc_ip{1},Pbase,0);
        uexs = load([Rpar.rootn 'data\phaseStimSave\stim_tmp_' sprintf('%3.f',1000*Rpar.IntP.phaseStim.phaseshift)],'uexs');
        pU{p} = uexs.uexs(Rpar.IntP.phaseStim.sensStm(2),round(Rpar.obs.brn*(1/R.IntP.dt))+1:end);
        disp([EpsN p])
    end
    rmdir([R.rootn 'data\phaseStimSave\'],'s')
    % Create save version
    feat_sim_save{1,EpsN} = feat_sim_base; feat_sim_save{2,EpsN} = feat_sim_stim;
    xsim_ip{1,EpsN} = xsim_ip_base; xsim_ip{2,EpsN} = xsim_ip_stim;
    pU_save{EpsN} = pU;

        ip = ip+1;

%             figure(100)
%             subplot(2,2,ip)
%             plot(R.frqz,squeeze(feat_sim_save{1,1}{1}(1,4,4,1,:)))
%             hold on
%             plot(R.frqz,squeeze(feat_sim_save{2,1}{6}(1,4,4,1,:)))
%             set(gca, 'YScale', 'log');
%             figure(200)
%             subplot(4,1,ip)
%             plot( pU_save{state}{1}' )
%
%             figure(300)
%             subplot(4,1,ip)
%             plot([0 Rout.IntP.tvec_obs],xsim_ip{1,state}{1}{1}(4,:)' )
%             hold on
%             plot([0 Rout.IntP.tvec_obs],xsim_ip{2,state}{6}{1}(4,:)' )
%

    rootan = [Rorg.rootn 'data\phaseLockedStimEpsComp'];
    mkdir(rootan)

    save([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_feat' num2str(SScomb) '.mat'],'feat_sim_save')
    save([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_xsim' num2str(SScomb) '.mat'],'xsim_ip','-v7.3')
    save([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_Rout' num2str(SScomb) '.mat'],'Rout')
    save([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_' num2str(EpsN) '_pU_save' num2str(SScomb) '.mat'],'pU_save')
end


% ! shutdown /s
