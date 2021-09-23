function fingerprintCompareWithEmpirical(R)
scmap = [0 0 0; brewermap(4,'Set1'); brewermap(2,'Set2')];
close all
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

% Setup State
statelistSpont = [          % state CON
    1 1;    % Fitted
    2 1;    % HD Down
    3 1;    % HD UP
    2 2;    % PS DOWN
    3 2;    % PS UP
    4 1;    % Control data    
    5 1];   % Lesion Data
% connectMat{band,NetState,CON}
cmat = []; smat = [];
for state = 1:size(statelistSpont,1)
    cmat(:,:,state) = [abs(connectMat{1,statelistSpont(state,1),statelistSpont(state,2)}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2)}),...
        angle(connectMat{1,statelistSpont(state,1),statelistSpont(state,2)}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2)})];
    smat(:,:,state) = specMat{1,statelistSpont(state,1),statelistSpont(state,2)};
    omat(:,:,state) = [statBurstOverl{1,statelistSpont(state,1),statelistSpont(state,2)}(:,:,1)+statBurstOverl{2,statelistSpont(state,1),statelistSpont(state,2)}(:,:,1)];
end


%% Stim Matrices
% Now Load in Stim Data
rootan = [R.rootn 'data\phaseLockedStim'];
load([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysis_CON_' num2str(1) '_feat' num2str(1) '.mat'],'connectMat','specMat','statBurstOverl')

% connectMat{band,phi,stm}
statelistStim = 1:13;
cmatStim = []; smatStim = [];
for state = 1:size(statelistStim,2)
    if state == 1
        phi = 1; stm = 1;
    else
        phi = state-1; stm = 2;
    end
    cmatStim(:,:,state) = [abs(connectMat{1,phi,stm}+connectMat{2,phi,stm}) angle(connectMat{1,phi,stm}+connectMat{2,phi,stm})];
    smatStim(:,:,state) = specMat{1,phi,stm};
    omatStim(:,:,state) = [statBurstOverl{1,phi,stm}(:,:,1)+statBurstOverl{2,phi,stm}(:,:,1)];
end

%% Do comparison
nmlist = [];
for stateStim = 1:size(cmatStim,3)
    A{2} = squeeze(cmatStim(:,:,stateStim));
    A{1} = normalize(squeeze(smatStim(:,:,stateStim)),'zscore','std');
    %     A{3} =  squeeze(omatStim(:,:,stateStim));
    for stateSpont = 1:size(cmat,3)
        B{2} = squeeze(cmat(:,:,stateSpont));
        B{1} =  normalize(squeeze(smat(:,:,stateSpont)),'zscore','std');
        
        %         B{3} = squeeze(omat(:,:,stateSpont));
        [h pnm(stateStim,stateSpont)] = corr(A{2}(:),B{2}(:),'type','Spearman');
        Lpart{1} = A{2}(1:6,1:6); Cpart{1} = A{2}(1:6,7:12); % seperate out linear vs circular parts of A
        Lpart{2} = B{2}(1:6,1:6); Cpart{2} = B{2}(1:6,7:12); % seperate out linear vs circular parts of B
        [Cpart{1},Cpart{2}] = remnan(Cpart{1}(:),Cpart{2}(:));
%         nmlist(stateStim,stateSpont) =  rsquare(A{2}(:),B{2}(:));
        nmlist(stateStim,stateSpont) = mean([rsquare(Lpart{1}(:),Lpart{2}(:)) circ_corrcc(Cpart{1}(:),Cpart{2}(:)).^2]); %rsquare(A{2}(:),B{2}(:));%*(pnm(stateStim,stateSpont)<0.05) ;% norm(A-B,'fro'); %%
        [h p] = corr(A{1}(:),B{1}(:),'type','Spearman');
        slist(stateStim,stateSpont) =  rsquare(A{1}(:),B{1}(:)); %*(p<0.05);%
        %         [h p] = corr(A{3}(:),B{3}(:),'type','Spearman');
        %         olist(stateStim,stateSpont) =  rsquare(A{3}(:),B{3}(:));%*(p<0.05);%
        clist(stateStim,stateSpont) =  collectR2(A,B);
        %         olist(st
    end
    
end



figure
plotFingerPrintMatch(nmlist,slist,clist,scmap,1:5,{'Fitted','HD-Down','HD-Up','PS-Down','PS-Up'},0)
subplot(1,3,1); ylim([-0.5 1]); subplot(1,3,2); ylim([0 1]); subplot(1,3,3); ylim([0 1]);

figure
plotFingerPrintMatch(nmlist,slist,clist,scmap,6:7,{'Cont.','Lesion'},1)
% subplot(1,3,1); ylim([-0.5 1]); subplot(1,3,2); ylim([0 1]); subplot(1,3,3); ylim([0 1]);

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
        %     lm = 100.*(lm-min(lm))./min(lm);
%         lm = 100*(lm - mean(lm));
%         ytit = {'Change in '; 'explained variance'};
        ytit = {'Similarity to '; 'Spontaneous (R2)'};
    else
        ytit = {'Similarity to '; 'Spontaneous (R2)'};
    end
    
    
    subplot(1,3,L)
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    phaseShiftDeg = rad2deg(phaseShift);
    p = plot(phaseShiftDeg,lm(2:13,netstates),'LineWidth',2);
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

