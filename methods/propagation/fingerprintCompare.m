function fingerprintCompare(R)
scmap = [0 0 0; brewermap(4,'Set1')];

%% Spontaneous Matrices
% Load Canonical States (predefined)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
load([rootan '\BB_' R.out.tag '_stateConnectivityMatrix.mat'],'connectMat','specMat','statBurstOverl')

% Setup State
statelistSpont = [1 1; 2 1; 3 1; 2 3; 3 3]; % state CON
% connectMat{band,state,CON}
cmat = []; smat = [];
for state = 1:size(statelistSpont,1)
    cmat(:,:,state) = [abs(connectMat{1,statelistSpont(state,1),statelistSpont(state,2)}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2)}),...
        angle(connectMat{1,statelistSpont(state,1),statelistSpont(state,2)}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2)})];
    smat(:,:,state) = specMat{1,statelistSpont(state,1),statelistSpont(state,2)};
    omat(:,:,state) = [statBurstOverl{1,statelistSpont(state,1),statelistSpont(state,2)}(:,:,1)+statBurstOverl{2,statelistSpont(state,1),statelistSpont(state,2)}(:,:,1)];
end

% Now Load in Stim Data
rootan = [R.rootn 'data\' R.out.oldtag '\phaseLockedStim'];
load([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysis_CON_' num2str(1) '_feat' num2str(1) '.mat'],'connectMat','specMat','statBurstOverl')

%% Stim Matrices
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
statelistStim = 1:13;
for stateStim = 1:size(cmatStim,3)
    A{2} = squeeze(cmatStim(:,:,stateStim));
    A{1} = normalize(squeeze(smatStim(:,:,stateStim)),'zscore','std');
%     A{3} =  squeeze(omatStim(:,:,stateStim));
    for stateSpont = 1:size(cmat,3)
        B{2} = squeeze(cmat(:,:,stateSpont));
        B{1} =  normalize(squeeze(smat(:,:,stateSpont)),'zscore','std');;
%         B{3} = squeeze(omat(:,:,stateSpont));
        [h pnm(stateStim,stateSpont)] = corr(A{2}(:),B{2}(:),'type','Spearman');
        nmlist(stateStim,stateSpont) =  rsquare(A{2}(:),B{2}(:))*(pnm(stateStim,stateSpont)<0.05) ;% norm(A-B,'fro'); %%
        [h p] = corr(A{1}(:),B{1}(:),'type','Spearman');
        slist(stateStim,stateSpont) =  rsquare(A{1}(:),B{1}(:)); %*(p<0.05);%
%         [h p] = corr(A{3}(:),B{3}(:),'type','Spearman');
%         olist(stateStim,stateSpont) =  rsquare(A{3}(:),B{3}(:));%*(p<0.05);%
        clist(stateStim,stateSpont) =  collectR2(A,B);
        %         olist(st
    end
   
end

figure
for L = 1:4
    if L == 2
        lm = nmlist;
    elseif  L == 1
        lm = slist;
    elseif L == 3
        lm = olist;
    elseif L == 4
        lm = clist;
    end
    
    subplot(1,4,L)
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    phaseShiftDeg = rad2deg(phaseShift);
    p = plot(phaseShiftDeg,lm(2:13,:),'LineWidth',2);
    for pip = 1:numel(p)
        p(pip).Color = scmap(pip,:);
    end
    xlabel('Stimulation Phase (radians)');
    xlim([0 360])
    grid on; box off
end

subplot(1,3,1)
title('Mapping to State Spectra')
ylabel('Similarity to Spontaneous (R2)')
ylim([-1 1])

subplot(1,4,2)
title('Mapping to State PLV')
ylabel('Similarity to Spontaneous (R2)')
ylim([0 1])

subplot(1,4,3)
title('Mapping to State Burst Overlap')
ylabel('Similarity to Spontaneous (R2)')
ylim([-1 1])
subplot(2,4,4)
title('Overall Mapping')
ylabel('Similarity to Spontaneous (R2)')
legend({'Base','HD-Down','HD-Up','PS-Down','PS-Up'},'Location','southeast');
ylim([-1 1])

