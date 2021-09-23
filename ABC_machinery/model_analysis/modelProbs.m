function [permMod, xsimMod] = modelProbs(x,m,p,R)
if ~isfield(R.analysis.BAA,'flag')
    R.analysis.BAA.flag = 0;
end
% load([R.rootn 'outputs\' R.out.tag '\parBank_' R.out.tag '_' d '.mat'])
parOptBank = R.parOptBank;

% parOptBank = parOptBank(:,parOptBank(end,:)>eps);
%% Compute KL Divergence
[KL DKL] = KLDiv(R,p,m,1);
R.Mfit.DKL = DKL;
N = R.analysis.modEvi.N;


    base = parOptBank(1:end,:);
    % If doing BAA analysis of model
    parforArg =0;
    ppm = [];
    N = 1;
    
    switch R.analysis.BAA.redmeth
        case 'average'
            % Take the expected parameters from distribution
            par = [];
            par{1} = spm_unvec(mean(base,2),p);
        case 'best'
            par = [];
            par{1} = spm_unvec(base(:,end),p);
        case 'UQ'
            %             nmrse = base(end-1,:);
            nmrse = base(end,:);
            X = base(:,nmrse>prctile(nmrse,75));
            par = [];
            par{1} = spm_unvec(mean(X,2),p);
    end

%%
figure(5)
pnew = par{1};
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
[r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData120319(R,m,u,pnew,0,1);
wfstr = ones(1,N);
while wfstr(end)>0
    parfor (jj = 1:N, parforArg)
        %     parfor jj = 1:N
        pnew = par{jj};
        u = innovate_timeseries(R,m);
        u{1} = u{1}.*sqrt(R.IntP.dt);
        [r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData120319(R,m,u,pnew,0);
        [ACC R2w] = computeObjective(R,r2)
        %     R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
        wfstr(jj) = any(wflag);
        r2rep{jj} = r2;
        accrep{jj} = ACC;
        par_rep{jj} = pnew;
        feat_rep{jj} = feat_sim;
        disp(jj); %
        if ~R.analysis.BAA.flag
            ppm.increment();
            xsims_rep{jj} = [];
        else
            xsims_rep{jj} = xsims_gl;
        end
    end
    
    if ~R.analysis.BAA.flag
        wfstr(end) = 0;
    end
end
delete(ppm);
permMod.r2rep = r2rep;
permMod.par_rep = par_rep;
permMod.feat_rep = feat_rep;
permMod.DKL = DKL;
permMod.KL = KL;
permMod.ACCrep = accrep;
xsimMod = xsims_rep;
permMod.MAP = spm_unvec(median(base,2),p);
[a b] = max([r2rep{:}]);
permMod.bestP = spm_unvec(base(:,b),p);

