function fingerprintCompare_Sweep_MSEALT(R)
close all
% GET EMPIRICAL TO PLOT AT SAME TIME!
scmap = [0 0 0; brewermap(4,'Set1'); brewermap(2,'Set2')];
% scmap = [brewermap(4,'Set1')];
% scmap = linspecer(5);
%% Spontaneous Matrices
% Load Canonical States (predefined)
rootan = [R.rootn 'data\ConnectionSweep'];
load([rootan '\BB_' R.out.tag '_stateConnectivityMatrix.mat'],'connectMat','specMat','statBurstOverl')
load([rootan '\BB_' R.out.tag '_EmpiricalStates_ConnectivityMatrix.mat'],'specMatEmp','connectMatEmp','statBurstOverlEmp')
% Add Empirical states (overwrite state 4 which is a permutation)
% You need to block out GPi/Thal comparisons for empirical states
connectMat(1:2,4:5,1) = connectMatEmp;
specMat(1:2,4:5,1) = specMatEmp;
statBurstOverl(1:2,4:5,1) = statBurstOverlEmp;
clear connectMatEmp specMatEmp statBurstOverlEmp

%% Get Connection Strengths for each CON
ck_1(1,:) = [1 logspace(-1,log10(5),34)];
ck_1(2,:) = [1 logspace(-1,log10(1.90),34)];

% Setup State
statelistSpont = [          % state CON
    1 1;    % Fitted
    2 1;    % HD Down
    3 1;    % HD UP
    2 2;    % PS DOWN
    3 2;    % PS UP
    4 1;    % Control data
    5 1];   % Lesion Data

complist = 2:5;
complistF = 2:7;
% connectMat{band,state,CON}
stateDefPLV = []; smat = [];
for state = 1:size(statelistSpont,1)
    stateDefPLV(:,:,state) = [abs(connectMat{1,statelistSpont(state,1),statelistSpont(state,2)}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2)}),...
        angle(connectMat{1,statelistSpont(state,1),statelistSpont(state,2)}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2)})];
    stateDefSpec(:,:,state) = specMat{1,statelistSpont(state,1),statelistSpont(state,2)};
    stateDefbOv(:,:,state) = [statBurstOverl{1,statelistSpont(state,1),statelistSpont(state,2)}(:,:,1)+statBurstOverl{2,statelistSpont(state,1),statelistSpont(state,2)}(:,:,1)];
end
%% Stim data

CON = 2; SScomb = 1; stm = 2;
for state = 1:31
    % Now Load in Stim Data
    rootan = [R.rootn 'data\phaseLockedStim'];
    load([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(state) '_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
        'connectMat','specMat','statBurstOverMat');
    %      load([rootan '\BB_InDrt_ModCompRev2_phaseLockedStim_burstAnalysis_CON_1_feat1.mat'],...
    %         'connectMat','specMat','statBurstOverl');
    statBurstOverMat = statBurstOverl;
    %% Stim Matrices
    % connectMat{band,phi,stm}
    philistStim = 1:12;
    %     cmatStim = []; smatStim = [];
    for phi = 1:size(philistStim,2)
        matPLVStim(:,:,phi,state) = [abs(connectMat{1,phi,stm}+connectMat{2,phi,stm}) angle(connectMat{1,phi,stm}+connectMat{2,phi,stm})];
        matSPECStim(:,:,phi,state) = specMat{1,phi,stm};
        %         matBOvStim(:,:,phi,state) = [statBurstOverMat{1,phi,stm}(:,:,1)+statBurstOverMat{2,phi,stm}(:,:,1)];
    end
end
%% Do comparison
nmlist = []; plvR2 = []; specR2 = [];
for state = 1:31
    for phiStim = 1:size(philistStim,2)
        A{2} = squeeze(matPLVStim(:,:,phiStim,state));
        A{1} =  normalize(squeeze(matSPECStim(:,:,phiStim,state)));
        subplot(3,1,1); plot(A{1}(:)); hold on
        %         A{3} =  squeeze(matBOvStim(:,:,phiStim,state));
        for stateSpont = 1:size(stateDefPLV,3)
            B{2} = squeeze(stateDefPLV(:,:,stateSpont));
            B{1} =  normalize(squeeze(stateDefSpec(:,:,stateSpont)));
            
            % Plot Stimulated
            subplot(3,1,2);  plot(B{1}(:)); hold on
            
            % Seperate linear and circular parts
            Lpart{1} = A{2}(1:6,1:6); Cpart{1} = A{2}(1:6,7:12); % seperate out linear vs circular parts of A
            Lpart{2} = B{2}(1:6,1:6); Cpart{2} = B{2}(1:6,7:12); % seperate out linear vs circular parts of B
            [Cpart{1},Cpart{2}] = remnan(Cpart{1}(:),Cpart{2}(:));
            [Lpart{1},Lpart{2}] = remnan(Lpart{1}(:),Lpart{2}(:));
            plvR2(stateSpont) = mean([-fxRMSE(Lpart{1}(:),Lpart{2}(:)) -fxCRMSE(Cpart{1}(:),Cpart{2}(:))]);
            
            specR2(stateSpont) =  -fxRMSE(A{1}(:),B{1}(:)); subplot(3,1,3);
            
            % Plot concat spectra
            subplot(3,1,3);  plot(A{1}(:)-B{1}(:)); hold on
            combR2(stateSpont) =  sum([ specR2(stateSpont)  plvR2(stateSpont)]);
        end
        plvStateList(state,phiStim) = complist(plvR2(complist)==max(plvR2(complist)));
        specStateList(state,phiStim) = complist(specR2(complist)==max(specR2(complist)));
        %         bOvStateList(state,phiStim) = find(bOvR2==max(bOvR2));
        combStateList(state,phiStim) = complist(combR2(complist)==max(combR2(complist)));
        
        plvR2Store(:,state,phiStim) = (plvR2(complistF));
        specR2Store(:,state,phiStim) = (specR2(complistF));
        %         bOvR2Store(:,state,phiStim) = (bOvR2);
        combR2Store(:,state,phiStim) = (combR2(complistF));
    end
end

figure
plotFingerPrintMatch(squeeze(plvR2Store(:,1,:)),squeeze(specR2Store(:,1,:)),...
    squeeze(combR2Store(:,1,:)),scmap,2:5,{'HD-Down','HD-Up','PS-Down','PS-Up'},0)
% subplot(1,3,1); ylim([-0.5 1]); subplot(1,3,2); ylim([-0.5 1]); subplot(1,3,3); ylim([0 1]);

figure

plotFingerPrintMatch(squeeze(plvR2Store(:,1,:)),squeeze(specR2Store(:,1,:)),...
    squeeze(combR2Store(:,1,:)),scmap,5:6,{'Cont.','Lesion'},1);

for i = 1:3
    switch i
        case 1 %
            AB = plvStateList;
            AB_r2 = plvR2Store;
        case 2
            AB = specStateList;
            AB_r2 = specR2Store;
            %         case 3
            %             AB = bOvStateList;
            %             AB_r2 = bOvR2Store;
        case 3; %4
            AB = combStateList;
            AB_r2 = combR2Store;
    end
    figure(100)
    subplot(2,3,i+3)
    
    [lr,zeroind] = min(abs(ck_1(CON,2:31)-1));
    if (ck_1(CON,zeroind+1)-1)<0
        list = [2:zeroind-1 1 zeroind:31];
    else
        list = [2:zeroind 1 zeroind+1:31];
    end
    %             zeroind = zeroind -1;
    
    plotMats(ck_1(CON,list),AB(list,1:end),scmap(complist,:),zeroind)
    
    
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    phaseShiftDeg = rad2deg(phaseShift);
    p = plot(phaseShiftDeg,squeeze(AB_r2(:,1,:)),'LineWidth',2); hold on
    for pip = 1:numel(p)
        p(pip).Color = scmap(complist(pip),:);
    end
    xlabel('Stimulation Phase (radians)');
    xlim([0 330]); %ylim([-0.1 1])
    grid on; box off;  axis square;
    a = gca;
    a.XTick = phaseShiftDeg(1:3:12);
    a.XTickLabel =  rad2deg(phaseShift(1:3:12));
    
end
set(gcf,'Position',[323          95        1435         883])

a = 1;
function LEG = plotFingerPrintMatch(nmlist,slist,clist,scmap,netstates,legnames,percflag)
for L = 1:3
    if L == 2
        lm = nmlist;
    elseif  L == 1
        lm = slist;
        %     elseif L == 3
        %         lm = olist;
    elseif L == 3
        lm = clist;
    end
    
    if percflag ==1
        %             lm = 100.*(lm-min(lm))./min(lm);
        lm = 100*(lm - mean(lm,2));
        ytit = {'Change in '; 'explained variance'};
        %         ytit = {'Similarity to '; 'Spontaneous (R2)'};
    else
        ytit = {'Similarity to '; 'Spontaneous (R2)'};
    end
    
    
    subplot(1,3,L)
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    phaseShiftDeg = rad2deg(phaseShift);
    p = plot(phaseShiftDeg,lm(netstates,:),'LineWidth',2);
    for pip = 1:numel(p)
        p(pip).Color = scmap(netstates(pip),:);
    end
    xlabel('Stimulation Phase (radians)');
    xlim([0 360])
    grid on; box off; axis square
end

subplot(1,3,1)
title('Mapping to State Spectra')
ylabel(ytit)

subplot(1,3,2)
title('Mapping to State PLV')
ylabel(ytit)

% subplot(1,4,3)
% title('Mapping to State Burst Overlap')
% ylabel('Similarity to Spontaneous (R2)')
% ylim([-1 1])

subplot(1,3,3)
title('Overall Mapping')
ylabel(ytit)
set(gcf,'Position',[326 239 1091 523])

LEG = legend(legnames);
LEG.Orientation = 'horizontal';
LEG.Position = [0.304 0.0634 0.437 0.0568];
LEG.Box = 'off';
LEG.Color = 'none';


function a = plotMats(y,AB,cmap,zeroind)
phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12

imagesc(AB,'AlphaData',~isnan(AB))
colormap(gca,cmap)
a = gca;
set(a, 'ydir', 'normal');
% c = colorbar;
% ylabel(c, 'Best Matching State')

caxis([2 5])
a.XTick = 1:3:12;
a.XTickLabel =  rad2deg(phaseShift(1:3:12));
a.YTick = [fliplr(zeroind:-4:1) zeroind+4:size(y,2)];
L = y([fliplr(zeroind:-4:1) zeroind+4:size(y,2)])*100;
a.YTickLabel = arrayfun(@(a) mat2str(a,3),L,'UniformOutput',0);
% a.XLabel.String = 'Sensing';
% a.YTick = 1:6;
% a.YTickLabel =     R.chsim_name;
% a.YLabel.String = 'Overlapping';
% hold on
% plot([7 0],[7 0],'k','LineWidth',2)
% c.Ticks = 1:4
axis square;
grid on



