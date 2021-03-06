function [R] = BAA_sim_phaseLockedStim(Rorg)

% Comopute simulations by sweeping across data
% [R,m,permMod] = getSimModelData_v3(Rorg,modID,simtime);
% mkdir([Rorg.rootn 'data\ModelFit\'])
% save([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')

% OR Load it in:
load([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')
for SScomb = 3; %1:2
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
    elseif SScomb == 3
        % Stimulating  STN
        senssite = 3; % GPe
        stimsite = 1; % STN
        stim_sens = 'stimSTN_sensGPe';        
    end
    R.IntP.phaseStim.sensStm = [senssite stimsite];
    
    %% Connection Sets
    ck_1(1,:) = [1 logspace(-1,log10(5),34)];
    ck_1(3,:) = [1 logspace(-1,log10(1.90),34)];
    
%     ck_1(1,:) = [1 ck_1_org(1,[2 31])];
%     ck_1(3,:) = [1 ck_1_org(3,[2 31])];
    
%     % Simulation Coniditions
    R.obs.csd.df = 0.5;
    R = setSimTime(R,128);
    
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
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    
    
    %% Loop through Connections
    for CON = [1 3]
        feat_sim_save = {};
        xsim_ip = {};
        for state = 1:size(ck_1,2)
            %% Setup Base Model
            Pbase = XBase;
            if CON == 1 % Hyperdirect
                Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(CON,state)); %
            elseif CON == 3 % Pallidal-subthalamo
                Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(CON,state)); %
            end
            
            % Setup stim parameterss
            R = typeIstimPars_v3(R);
            
            % Simulate Base Model
            uc_ip{1} = uc;
            R.IntP.phaseStim.switch = 0 ;
            R.IntP.phaseStim.phaseshift = 0;
            [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1},~,Rout] = computeSimData(R,m,uc_ip{1},Pbase,0);
            
            % Work out the threshold
           [~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);
            eps(state) =  R.IntP.phaseStim.eps;
            
            %% Do Stimulation
            R.IntP.phaseStim.switch = 1;
            m = m; % initialise for parfor
            xsim_ip_stim = cell(1,12); feat_sim_stim = cell(1,12); pU = cell(1,12);
            parfor p = 1:numel(phaseShift)
                Rpar = R;
                % Modulate the phase
                Rpar.IntP.phaseStim.phaseshift = phaseShift(p);
                
                % Simulate with Stimulation
                [~,~,feat_sim_stim{p},~,xsim_ip_stim{p},~,Rpar]  = computeSimData(Rpar,m,uc_ip{1},Pbase,0);
                uexs = load([Rpar.rootn 'data\rat_InDirect_ModelComp\phaseStimSave\stim_tmp_' sprintf('%3.f',1000*Rpar.IntP.phaseStim.phaseshift)],'uexs');
                pU{p} = uexs.uexs(Rpar.IntP.phaseStim.sensStm(2),round(Rpar.obs.brn*(1/R.IntP.dt))+1:end);
                disp([CON state p])
            end
            rmdir([R.rootn 'data\rat_InDirect_ModelComp\phaseStimSave\'],'s')
            % Create save version
            feat_sim_save{1,state} = feat_sim_base; feat_sim_save{2,state} = feat_sim_stim;
            xsim_ip{1,state} = xsim_ip_base; xsim_ip{2,state} = xsim_ip_stim;
            pU_save{state} = pU;
        end
        rootan = [Rorg.rootn 'data\' Rorg.out.oldtag '\phaseLockedStim'];
        mkdir(rootan)
        
        save([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],'feat_sim_save')
        save([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip','-v7.3')
        save([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_Rout' num2str(SScomb) '.mat'],'Rout')
        save([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_pU_save' num2str(SScomb) '.mat'],'pU_save')
    end
    
end
% ! shutdown /s