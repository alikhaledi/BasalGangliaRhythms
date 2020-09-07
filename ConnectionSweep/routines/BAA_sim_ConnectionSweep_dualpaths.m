function [R] = BAA_sim_ConnectionSweep_dualpaths(Rorg,modID,simtime,HD)

% Comopute simulations by sweeping across data
% [R,m,permMod] = getSimModelData_v3(Rorg,modID,simtime);
% mkdir([Rorg.rootn 'data\ModelFit\'])
% save([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')

% OR Load it in:
load([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
warning('Loading Preloaded model, cant change simtime or model choice!!!')
pause(1)
%% Connection Sets
if HD == 1
    % discrete connection list
    for CON = [1 3]
        ck_1(CON,:) = [0.00001 0.125 0.25 0.5 0.75 1 1.25  1.5 3 5];
    end
    hdext = '';
elseif HD == 2
    % Sliding scale
    % Used for: (1) Plot Spectra over sweeps
    for CON = [1 3]
        ck_1(CON,:) = logspace(-1,1,40);
    end
    hdext = '_F1';
elseif HD == 3
    % This compute the ranges already used with the lower definition
    % simulation sweep in plotSweepSpectra.betaKrange(:,CON) = [b1 b2 b3];
    % i.e. [10% 100% 190%]. This used HD == 2 to setup the vector. Now
    % truncate using simulated beta values
    % This is currently used for burst statistics
    rootan = [Rorg.rootn 'data\rat_InDirect_ModelComp\ConnectionSweep'];
    load([rootan '\BB_' Rorg.out.tag '_ConnectionSweep_CON_KRange.mat'],'Krange')
    
    
%     refK = logspace(-1,0.7,30);
    for CON = [1 3]
        linspace(exp(Krange{CON}(1)),exp(Krange{CON}(end)),20)
            ck_1(CON,:) = linspace((Krange{CON}(1)),(Krange{CON}(end)),25);
            [a b] = min(abs(ck_1(CON,:)-1))
             ck_1(CON,b) = 1; % Ensure base model is included
%         x1 = logspace(log10(refK(betaKrange(1,CON))),log10(refK(betaKrange(2,CON))),10);
%         x2 = logspace(log10(refK(betaKrange(2,CON))),log10(refK(betaKrange(3,CON))),10); 
%         ck_1(CON,:) = [x1(1:end-1) 1 x2(2:end)];
    end
    hdext = '_bKF';
    
    elseif HD == 4
    % Sliding scale, deterministic
    % Used for: (1) Plot Spectra over sweeps
    for CON = [1 3]
        ck_1(CON,:) = logspace(-1,0.7,15);
    end
    hdext = '_F1';
    m.uset.p.scale = 0;
    
end

% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
XBase = permMod{1}.par_rep{1};
R.obs.SimOrd = 10;
R.obs.trans.norm = 0;
R.obs.gainmeth = {};

%% Loop through Connections
R.IntP.getNoise = 0;
m = m;
for CON = 3
    feat = {};
    xsim = {};
    for i = 1:size(ck_1,2)
        % Now Modify
        Pbase = XBase;
        if CON == 1 % Hyperdirect
            Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(CON,i)); %
            Pbase.A{2}(4,3) = -32; %

        elseif CON == 3 % Pallidal-subthalamo
            Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(CON,i)); %
            Pbase.A{1}(4,1) = -32; %
        end
        [r2mean,pnew,feat_sim,dum1,xsim_gl] = computeSimData(R,m,uc,Pbase,0);
        feat{i} = feat_sim;
        xsim{i} = xsim_gl;
    end
    rootan = [Rorg.rootn 'data\' Rorg.out.oldtag '\ConnectionSweep'];
    mkdir(rootan)
    
    save([rootan '\BB_' Rorg.out.tag '_ConnectionSweepDualPaths_CON_' num2str(CON) '_feat' hdext '.mat'],'feat')
    save([rootan '\BB_' Rorg.out.tag '_ConnectionSweepDualPaths_CON_' num2str(CON) '_xsim' hdext '.mat'],'xsim')
    save([rootan '\BB_' Rorg.out.tag '_ConnectionSweepDualPaths_CON_' num2str(CON) '_ck_1' hdext '.mat'],'ck_1')
end


