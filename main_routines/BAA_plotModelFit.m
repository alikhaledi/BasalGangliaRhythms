function [R] = BAA_plotModelFit(Rorg,modID,simtime)
% Simulate the Base Model
Rorg.obs.trans.norm = 1; % Normalize the output spectra
Rorg.obs.trans.gauss = 1; % Smooth with sum of 3 gaussians
Rorg.obs.obsstates = [1:6]; % All 6 nodes are observed
Rorg.chloc_name = Rorg.chsim_name; % Ensure sim names match to output names

%% Call the simulator
[R,m,permMod,xsimMod] = getSimModelData_v3(Rorg,modID,simtime,1);
% Save for later use (avoids recalling simulator)
mkdir([Rorg.rootn 'data\modelfit'])
save([Rorg.rootn 'data\modelfit\SimModelData_M10.mat'],'R','m','permMod')

close all;
X = xsimMod{1}{1}{1};
for i = 1:6
    figure(1)
    plot(R.IntP.tvec_obs,X(i,:))
    hold on
    %     figure(2+i)
    %  pspectrum(X(i,:),1/R.IntP.dt,'spectrogram','TimeResolution',0.5,...
    % 'OverlapPercent',99,'MinTHreshold',-20);
    % ylim([0 0.098])
    
end
xlim([30 40])
set(gcf,'Position',[50         550        1686         436])
R = prepareRatData_NoGauss_Group_NPD(Rorg,0,0);

figure
plotABCSpectraOnly(R.data.feat_xscale,R.data.feat_emp,permMod{1}.feat_rep{1})
figure
npdplotter_110717({R.data.feat_emp},{permMod{1}.feat_rep{1}},R.data.feat_xscale,R,[],[])
for i = 1:6
    for j = 1:6
        if i~=j
            subplot(6,6,sub2ind([6 6],j,i))
            ylim([0 0.6]);
        end
    end
end

