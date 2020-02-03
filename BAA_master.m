%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
R = basalGangliaRhythms_AddPaths();

%% Setup
%% Set Routine Pars
R = setupBasalGangliaModel(R);
modID = 10;
%% Plot Model Fit
simlength = 128;
BAA_plotModelFit(R,modID,simlength);

% Lesion analysis
% BAA_sim_lesionExp_cortical_rhythms(R,modID,32,1)
simlength = 32;
fresh = 1;
BAA_sim_lesionExp(R,modID,simlength,fresh)

%% Plot Model Sweep Spectra
BAA_sim_ConnectionSweep_v2(R,modID,100,2)
R = plotSweepSpectraWrapper(R); % This doesnt just plot but also works out the vectors for spread
plotSweepSpectraWrapper_M2_SI(R); % Plot Cortical Spectra

%% ISSUES IN THESE
% Make Clear which sets of data are being used.
simulateBurstData(R); % basically just gets'_bKF' data, scaled to beta 10,100,190%
R = computeBurstWrapper(R);
%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R); %,BB,[15 17 20]); % [15 17 20] for STN_GPE [1 6 8]
OnsetEvolveAnalysisMaster(R)

%% Compute Closed Loop Experiment
simlength = 50;
fresh = 0;
[R] = BAA_sim_fakeCloseLoop(R,modID,simlength,fresh);

% This is a analysis of Beta inputs STN/M2
simlength = 32;
fresh = 1;
[R,MP] = BAA_sim_betaInputs(R,modID,simlength,fresh);

% This is pertubation analysis to get PRCs- at the moment its not very
% clear what were looking for with it
% [R] = BAA_sim_PRC(R,MP,500,0);

[R] = BAA_sim_BetaPropagation(R,160,1); %remember reference of Boba and Hamacher 2015 for ZTE
 
 
% BB.range.RP = linspace(-pi,pi,7);
% BB = computeBetaBurstRPStats(R,BB);
% 
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)
[R,m,permMod,xsimMod] = getSimModelData_Draw(R,modID,32,0,150);
% This Draw is not working!

save('permModtmp_inflated','permMod')
load('permModtmp_inflated','permMod')
netA = [];
netAbsA = []; bpowr = []; bpowr_br = []; bcohr = [];
for i = 1:numel(permMod.wflag)
    if permMod.wflag(i)
        pinst = permMod.par_rep{i};
        feat = permMod.feat_rep{i};
        % X1) compute the "direct" loop - CTX->STR->GPi->Thal->CTX
        X1(1) = pinst.A{1}(2,1); % M2 -> STR
        X1(2) = -pinst.A{2}(5,2); % STR -> GPi
        X1(3) = -pinst.A{2}(6,5); % GPi -> Thal
        X1(4) = pinst.A{1}(1,6); % Thal -> CTX
        
        % X2) compute the "indirect" loop - CTX->STR->GPe->STN->GPi->Thal->CTX
        X2(1) = pinst.A{1}(2,1); % M2 -> STR
        X2(2) = -pinst.A{2}(3,2); % STR -> GPe
        X2(3) = -pinst.A{2}(4,3); % GPe -> STN
        X2(4) = pinst.A{1}(5,4); % STN -> GPi
        X2(5) = -pinst.A{2}(6,5); % GPi -> Thal
        X2(6) = pinst.A{1}(1,6); % Thal -> CTX
        
        % X3) compute the "hyperdirect" loop - CTX->STN->GPi->Thal->CTX
        X3(1) = pinst.A{1}(4,1); % M2 -> STN
        X3(2) = pinst.A{1}(5,4); % STN -> GPi
        X3(3) = -pinst.A{2}(6,5); % GPi -> Thal
        X3(4) = pinst.A{1}(1,6); % Thal -> CTX
        
        % X4) compute the "thalamocortical" loop - CTX->Thal->CTX
        X4(1) = pinst.A{1}(6,1); % M2 -> STN
        X4(2) = pinst.A{1}(1,6); % STN -> GPi
        
        % X5) compute the "pallido-subthalamic" loop - GPe->STN->GPe
        X5(1) = -pinst.A{2}(4,3); % GPe -> STN
        X5(2) = pinst.A{1}(3,4); % STN -> GPe
        
        % Compute Summaries
        netA(:,i) = [sum(X1) sum(X2) sum(X3) sum(X4) sum(X5)]';
        netAbsA(:,i) = [sum(abs(X1)) sum(abs(X2)) sum(abs(X3)) sum(abs(X4)) sum(abs(X5))]';
        
        
        [bpowr_br(i),fpow_br(i),bpowr(i),fpow(i),bcohr(i),fcoh(i)] = computeBetaSpectralStats(R.frqz,{feat})
    else
        netA(:,i) = nan(1,5);
        netAbsA(:,i) = nan(1,5);
        bpowr(i) = nan;
        bpowr_br(i) = nan;
        bcohr(i) = nan;
    end
end

for i = 1:5
    figure(i)
    scatterhist((netAbsA(i,:)),bcohr)
end

% COherence seems to give good correlations



for i = 1:100
plotABCSpectraOnly(R.data.feat_xscale,R.data.feat_emp,permMod.feat_rep{i})
end

for p = 1:4; subplot(1,6,p); ylim([0 1]); end

