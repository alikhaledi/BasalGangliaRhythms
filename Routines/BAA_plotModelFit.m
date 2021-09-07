function [R] = BAA_plotModelFit(Rorg,modID,simtime)
% Simulate the Base Model
% Rorg.obs.gainmeth = {'unitvar','boring'};
% Rorg.obs.trans.norm = 1; % Normalize the output spectra
% Rorg.obs.trans.gauss = 1; % Smooth with sum of 3 gaussians
% Call the simulator
 [Rout,m,p] = dataUpdateFix(Rorg,modID);
 Rout.obs.obsstates = [1:6]; % All 6 nodes are observed
Rout.chloc_name = Rout.chsim_name; % Ensure sim names match to output names
% Rorg.out.dag = sprintf([Rorg.out.tag '_M%.0f'],modID);
% [Rout,m,p] = loadABCData_160620(Rorg);
%  Rout = dataUpdateFix(Rorg,modID)
% Rout.IntP.intFx = @ABC_fx_compile_120319; % update fx name
[modelError,pnew,feat_sim,xsims,xsims_gl,~,Rout] = computeSimData_160620(Rout,m,[],p,simtime,1);
mkdir([Rorg.path.rootn 'data\modelfit'])
permMod{1}.par_rep{1} = p;%
save([Rorg.path.rootn 'data\modelfit\SimModelData_M10.mat'],'Rout','m','permMod')

close all;
X = xsims_gl{1};
for i = 1:6
    figure(1)
    plot(Rout.IntP.tvec_obs,X(i,2:end)-(6*i))
    hold on
    %     figure(2+i)
    %  pspectrum(X(i,:),1/R.IntP.dt,'spectrogram','TimeResolution',0.5,...
    % 'OverlapPercent',99,'MinTHreshold',-20);
    % ylim([0 0.098])
    
end
xlim([30 40])
set(gcf,'Position',[50         550        1686         436])
% R = prepareRatData_NoGauss_Group_NPD(Rorg,0,0);

figure
plotABCSpectraOnly(Rout.data.feat_xscale{1},Rout.data.feat_emp{1},feat_sim{1})
figure
genplotter_200420({Rout.data.feat_emp},{feat_sim},Rout.data.feat_xscale,Rout,[],[])

% featn = {'cross','auto','cross_only'};
% bandlim = {[],[14 21],[21 30],[14 30]};
% for band = 1:4
%     for feat = 1:3
%         R.objfx.bandlim = bandlim{band};
%         R.objfx.specspec = featn{feat};
%         
%         r2mean(band,feat)  = R.IntP.compFx(R,permMod{1}.feat_rep{1});
%     end
% end

%% SCRIPT GRAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:36
%     if ~any(intersect([1 8 15 22 29 36],i))
%         subplot(6,6,i)
%         ylim([0 1])
%     end
% end
% A = X(6,:);
%
% STN_slow = ft_preproc_bandpassfilter(X(4,:),2000,[0.5 4],[],'fir');
%
% slow_list = linspace(2.5,6,32);
% fast_list = linspace(8,34,32);
% Rz = [];
% for i = 1:size(slow_list,2)
%     for j = 1:size(fast_list,2)
%         A = ft_preproc_bandpassfilter(X(4,:),2000,[slow_list(i)-2 slow_list(i)+2],[],'fir');
%         B = abs(hilbert(ft_preproc_bandpassfilter(X(4,:),2000,[fast_list(j)-2 fast_list(j)+2],[],'fir')));
%
%        LR = corrcoef(A,B);
%        Rz(i,j) = LR(2);
%        if i==4
%            [r,lags] = xcorr(A,B,5000);
%            xcorstore(:,j) = r;
%        end
%         disp([i j])
%     end
% end
% subplot(1,3,1)
% imagesc(slow_list,fast_list,Rz')
% set(gca,'YDir','normal');
% subplot(1,2,2)
% imagesc(lags,fast_list,xcorstore')
% set(gca,'YDir','normal'); xlim([-1000 1000])
% % subplot(1,3,3)
% % imagesc(lags,fast_list,lagstore_j')
% % set(gca,'YDir','normal'); xlim([-1000 1000])

