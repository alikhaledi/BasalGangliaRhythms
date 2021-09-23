function fingerprintCompareSimulatedEmpirical(R)
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

%% Do comparison
nmlist = [];
empStateList = [6 7];
spontStateList = [1:5];
for empState = 1:size(empStateList,2)
    A{1} =  normalize(squeeze(smat(:,:,empStateList(empState))),'zscore','std');
    A{2} = squeeze(cmat(:,:,empStateList(empState))));
    for spontState = 1:size(spontStateList,2)
        B{1} =  normalize(squeeze(smat(:,:,spontStateList(spontState))),'zscore','std');
        B{2} = normalize(squeeze(cmat(:,:,spontStateList(spontState))));
        
        % Do the spectra first
        atmp = A{1}; btmp = B{1};
        [atmp, btmp] = remnan(atmp,btmp);
        [R p] = corr(atmp(:),btmp(:),'type','Spearman');
        Rsqr = rsquare(atmp(:),btmp(:));
                NRMSE = goodnessOfFit(atmp(:),btmp(:),'NRMSE'); 
        slist(:,empState,spontState) =  [Rsqr,R,NRMSE];
        
        % Now do the connectivity
        % Do the circular parts seperately
        Lpart{1} = A{2}(1:6,1:6); Cpart{1} = A{2}(1:6,7:12); % seperate out linear vs circular parts of A
        Lpart{2} = B{2}(1:6,1:6); Cpart{2} = B{2}(1:6,7:12); % seperate out linear vs circular parts of B
        [Cpart{1},Cpart{2}] = remnan(Cpart{1}(:),Cpart{2}(:));
                [Lpart{1},Lpart{2}] = remnan(Lpart{1}(:),Lpart{2}(:));
        Rsqr = mean([rsquare(Lpart{1}(:),Lpart{2}(:)) circ_corrcc(Cpart{1}(:),Cpart{2}(:)).^2]); %rsquare(A{2}(:),B{2}(:));%*(pnm(stateStim,stateSpont)<0.05) ;% norm(A-B,'fro'); %%
        R = mean([corr(Lpart{1}(:),Lpart{2}(:),'type','Spearman') circ_corrcc(Cpart{1}(:),Cpart{2}(:))]); %rsquare(A{2}(:),B{2}(:));%*(pnm(stateStim,stateSpont)<0.05) ;% norm(A-B,'fro'); %%
        NRMSE = goodnessOfFit(Lpart{1}(:),Lpart{2}(:),'NRMSE'); 
        nmlist(:,empState,spontState) =  [Rsqr,R,NRMSE];
         
        % Now do the complete set
        %         [h p] = corr(A{3}(:),B{3}(:),'type','Spearman');
        R = mean([corr(Lpart{1}(:),Lpart{2}(:),'type','Spearman') circ_corrcc(Cpart{1}(:),Cpart{2}(:)).^2]); %rsquare(A{2}(:),B{2}(:));%*(pnm(stateStim,stateSpont)<0.05) ;% norm(A-B,'fro'); %%
        clist(empState,spontState) =  collectR2(A,B);
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

