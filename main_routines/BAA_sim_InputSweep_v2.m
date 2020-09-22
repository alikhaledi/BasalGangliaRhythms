function [R] = BAA_sim_InputSweep_v2(Rorg,modID,simtime,HD)

% Comopute simulations by sweeping across data
% [R,m,permMod] = getSimModelData_v3(Rorg,modID,simtime);
% OR Load it in:
load([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')
warning('Loading Preloaded model, cant change simtime or model choice!!!')
pause(1)
    R = setSimTime(R,simtime);

XBase = permMod{1}.par_rep{1};
R.obs.SimOrd = 10;
R.obs.trans.norm = 0;
R.obs.gainmeth = {};

%% Loop through Noise Inputs
noiseVec = logspace(-3,1,30);
R.IntP.Utype = 'constant';
for np = 1:numel(noiseVec)
%     % Compute Noise
    m.uset.p.scale =noiseVec(np);
    uc = innovate_timeseries(R,m);
    uc{1} = uc{1}.*sqrt(R.IntP.dt);
    
    Pbase = XBase;
    [r2mean,pnew,feat_sim,dum1,xsim_gl] = computeSimData(R,m,uc,Pbase,0);
    nse = xsim_gl;
    
    % get process
    x = xsim_gl{1};
    % Find dx/dt ~= 0
    dstar = diff(x,1,2)<1e-32;
    for i = 1:size(x,1)
        dInds = find(dstar(i,:));
        if numel(dInds)>256
            dInds = dInds(randperm(numel(dInds),256));
        end
        dstarVal{i,np} = x(dstar(i,dInds));
%         dstarmean(i,np) = mean(x(i,dInds))
    end
    
    
end

for np = 1:numel(noiseVec)
    scatter(repmat(noiseVec(np),1,size(dstarVal{4,np},2)),dstarVal{4,np});
    hold on
end
