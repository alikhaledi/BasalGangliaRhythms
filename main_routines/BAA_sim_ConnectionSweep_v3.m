function [R] = BAA_sim_ConnectionSweep_v3(Rorg)

% MODID
load([Rorg.rootn 'data\ModelFit\SimModelData_M10.mat'],'R','m','permMod')
warning('Loading Preloaded model, cant change simtime or model choice!!!')
pause(1)

% Sliding scale
% Used for: (1) Plot Spectra over sweeps
ck_1(1,:) = [1 logspace(-1,log10(3.11),30)]; % This is HD
ck_1(2,:) = [1 logspace(-1,log10(1.33),30)]; % This is PS

hdext = '_REV';

%% Trans Options
R.obs.SimOrd = 10;
R.obs.trans.norm = 0;
R.obs.gainmeth = {};

%% Compute Noise
% Give all timeseries the same input - makes comparable
rng(2315324)
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
XBase = permMod{1}.par_rep{1};
R.IntP.getNoise = 1; % This turns off all the connections to get resting endogenous noise for each node
R.obs.csd.df = 0.25; % increase spectral resolution
R.obs.SimOrd = 11;
[dum1,dum2,feat_sim_noise,dum3,xsim_noise] = computeSimData(R,m,uc,XBase,0);
R.IntP.getNoise = 0;


%% Loop through Connections
m = m;
for CON = 1:2
    feat = {};
    xsim = {};
    parfor i = 1:size(ck_1,2)
        % Now Modify
        Pbase = XBase;
        if CON == 1 % Hyperdirect
            Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(CON,i)); %
        elseif CON == 2 % Pallidal-subthalamo
            Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(CON,i)); %
        end
        [r2mean,pnew,feat_sim,dum1,xsim_gl] = computeSimData(R,m,uc,Pbase,0,1);
        feat{i} = feat_sim;
        xsim{i} = xsim_gl;
        
        disp([CON i])
    end
    rootan = [Rorg.rootn 'data\' Rorg.out.tag '\ConnectionSweep'];
    mkdir(rootan)
    
    save([rootan '\BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext '.mat'],'feat')
    save([rootan '\BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim' hdext '.mat'],'xsim')
    save([rootan '\BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext '.mat'],'ck_1')
end

save([rootan '\BB_' Rorg.out.tag '_ConnectionSweep_noise.mat'],'xsim_noise','feat_sim_noise')

