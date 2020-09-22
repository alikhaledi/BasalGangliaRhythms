function [R] = BAA_sim_lesionExp_signChange(Rorg,modID,simtime,fresh)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v3(Rorg,modID,simtime);
p = permMod{1}.par_rep{1};
R = setSimTime(R,simtime); % This sets the simulatiomn time by modify R structure
R.Bcond = -1; % No contrasting conditions

% Plotting colors
cmap = brewermap(128,'RdBu');
ckeypow = linspace(-100,75,128); % These make a key for colorcoding the results
ckeyfrq = linspace(-10,10,128);

% Connection names
condname = {'Fitted','M2->STR','M2->STN','M2->Thal','STN->GPe','STR-|Gpe','STR->GPi','GPe-|STN','STN->GPi','GPi-|Thal','Thal->M2'};
% Connection spec [ex/inh src tar]
AIJ = {[],[1 2 1],[1 4 1],[1 6 1],[1 3 4],[2 3 2],[2 5 2],[2 4 3],[1 5 4],[2 6 5],[1 1 6]};

% Give all timeseries the same input
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
R.obs.trans.norm = 0; % No normalization of spectra

Pbase = p;
if fresh == 1
    for i = 1; %:size(condname,2)
        Pbase_i = Pbase; % Load parameter base
        if i~=1
            XL = Pbase_i.A{AIJ{i}(1)}(AIJ{i}(2),AIJ{i}(3));
            Pbase_i.A{AIJ{i}(1)}(AIJ{i}(2),AIJ{i}(3)) = -32; % Remove the ith connection
            if AIJ{i}(1) == 1
                Pbase_i.A{2}(AIJ{i}(2),AIJ{i}(3)) = log(exp(XL).*1); % Replace with Fliiped Sign
            elseif AIJ{i}(1) == 2
                Pbase_i.A{1}(AIJ{i}(2),AIJ{i}(3)) =log(exp(XL).*1); % Replace with Fliiped Sign
            end
        end
        parfor j = 1:size(condname,2)    % Setup the simulations
            Pbase_ij = Pbase_i;
            if j~=1
                YL = Pbase_ij.A{AIJ{j}(1)}(AIJ{j}(2),AIJ{j}(3));
                Pbase_ij.A{AIJ{j}(1)}(AIJ{j}(2),AIJ{j}(3)) = -32; % Remove the jth connection
                if AIJ{j}(1) == 1
                    Pbase_ij.A{2}(AIJ{j}(2),AIJ{j}(3)) = log(exp(YL).*1); % Flip Sign
                elseif AIJ{j}(1) == 2
                    Pbase_ij.A{1}(AIJ{j}(2),AIJ{j}(3)) = log(exp(YL).*1); % Flip Sign
                end
            end
            [r2,~,feat_sim] = computeSimData(R,m,uc,Pbase_ij,0); % Simulates the new model
            
            % Find resulting power statistics of the simulated data (STN
            % beta
            [powIJ_B(j,i),peakIJ_B(j,i),freqIJ_B(j,i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[14 30]);
            [powIJ_B1(j,i),peakIJ_B1(j,i),freqIJ_B1(j,i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[14 21]);
            [powIJ_B2(j,i),peakIJ_B2(j,i),freqIJ_B2(j,i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,1,:)),[21 30]);
            fitIJ(j,i) = r2;
            disp([i j])
        end
    end
    rootan = [Rorg.rootn 'data\' Rorg.out.oldtag '\BAA_InActivation'];
    mkdir(rootan)
    save([rootan '\BAA_InActivation'],'powIJ_B','peakIJ_B','freqIJ_B',...
        'powIJ_B1','peakIJ_B1','freqIJ_B1',...
        'powIJ_B2','peakIJ_B2','freqIJ_B2','condname')
else
    rootan = [Rorg.rootn 'data\' Rorg.out.oldtag '\BAA_InActivation'];
    load([rootan '\BAA_InActivation'],'powIJ_B','peakIJ_B','freqIJ_B',...
        'powIJ_B1','peakIJ_B1','freqIJ_B1',...
        'powIJ_B2','peakIJ_B2','freqIJ_B2','condname')
end
% Plot Results
figure
for band = 1
    titlist = {'14-30 Hz','14-21 Hz','21-30 Hz'}
    if band == 1 % Whole band 13-30HZ
        X = (powIJ_B-powIJ_B(1,1))./(powIJ_B(1,1)).*100;
        Y = (freqIJ_B-freqIJ_B(1,1));
    elseif band == 2 % Low beta band 14-21 Hz
        X = (powIJ_B1-powIJ_B1(1,1))./(powIJ_B1(1,1)).*100;
        Y = (freqIJ_B1-freqIJ_B1(1,1));
    elseif band == 3 % High beta band 21-30 Hz
        X = (powIJ_B2-powIJ_B2(1,1))./(powIJ_B2(1,1)).*100;
        Y = (freqIJ_B2-freqIJ_B2(1,1));
    end
    % Color Setup
    bcmap = brewermap(128,'RdGy');
    
    
    subplot(1,2,1)
    colormap(bcmap)
    igsc = imagesc(X)
    set(gca,'YDir','normal');
    cb = colorbar;
    caxis([-100 100])
    %         title(titn)
    
    % Plot Barplots for individual connections
    subplot(2,1,1)
    X = X';
    L = X(1,:);
    L(L>1e3) = NaN;
    b = bar(L);
    [dum keyind] = min(abs(L'-ckeypow),[],2);
    b.CData = cmap(keyind,:);
    b.EdgeAlpha = 0; b.FaceColor = 'flat';
    a = gca; grid on; box off
    a.XTickLabel = condname; a.XTickLabelRotation = 45; ylabel('% Power Change')
    xlim([1.5 11.5])
    
    subplot(2,1,2)
    Y = Y';
    L = Y(1,:);
    L(L>1e3) = NaN;
    b = bar(L);
    [dum keyind] = min(abs(L'-ckeyfrq),[],2);
    b.CData = cmap(keyind,:);
    b.EdgeAlpha = 0; b.FaceColor = 'flat';
    a = gca; grid on; box off
    a = gca;
    a.XTickLabel = condname; a.XTickLabelRotation = 45; ylabel('Change in Peak Frequency')
    xlim([1.5 11.5])
    
    % Make Tables
%     tabdata(:,band) = diag(X);
    
end
set(gcf,'Position',[448   121   559   592])

%% Heatmap of connection synergies
% Compute synergies
Z = (powIJ_B-powIJ_B(1,1))./(powIJ_B(1,1)).*100; % This is the simulate effect (IJ)
Z(getDiagInd(Z)) = zeros(1,size(Z,1));
X = powIJ_B;
Y = X;
for i = 1:size(X,1)
    for j= 1:size(X,2)
        if i~=j
            Y(i,j) = ((X(i,i)+X(j,j))-powIJ_B(1,1))./(powIJ_B(1,1)).*100; % This is the additive effect (I) + (J)
        end
    end
end

% plot heatmap
figure
colormap(flipud(cmap))
Q = abs(Y)-abs(Z)
Z(isnan(Z)) = 0;
igsc = imagesc(Q)
set(gca,'YDir','normal');
caxis([-125 125]);
a = gca; a.XTick = 1:size(X,1); a.XTickLabel = condname; a.XTickLabelRotation = 45;
a.YTickLabel = condname; a.YTickLabelRotation = 0;
cb = colorbar; axis equal; xlim([1.5 11.5]);ylim([0.5 11.5])

set(gcf,'Position',[1065          84         692         592])
