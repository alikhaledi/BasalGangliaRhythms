function [] = computeModelDrawParameters(R,permMod)


load('ConnectionStrengthEX.mat')

E(1) = EX(1,1); % MMC
E(2) = EX(2,2); % STR
E(3) = EX(3,2); % GPe
E(4) = EX(4,1); % STN
E(5) = EX(5,2); % GPi
E(6) = EX(6,1); % Thal.


netA = [];
netAbsA = []; bpowr = []; bpowr_br = []; bcohr = [];
for i = 1:numel(permMod.wflag)
    if permMod.wflag(i)
        pinst = permMod.par_rep{i};
        feat = permMod.feat_rep{i};
        % X1) compute the "direct" loop - CTX->STR->GPi->Thal->CTX
        X1(1) = E(1).*exp(pinst.A{1}(2,1)); % M2 -> STR
        X1(2) = E(2).*exp(pinst.A{2}(5,2)); % STR -> GPi
        X1(3) = E(5).*exp(pinst.A{2}(6,5)); % GPi -> Thal
        X1(4) = E(6).*exp(pinst.A{1}(1,6)); % Thal -> CTX
        
        % X2) compute the "indirect" loop - CTX->STR->GPe->STN->GPi->Thal->CTX
        X2(1) = E(1).*exp(pinst.A{1}(2,1)); % M2 -> STR
        X2(2) = E(2).*exp(pinst.A{2}(3,2)); % STR -> GPe
        X2(3) = E(3).*exp(pinst.A{2}(4,3)); % GPe -> STN
        X2(4) = E(4).*exp(pinst.A{1}(5,4)); % STN -> GPi
        X2(5) = E(5).*exp(pinst.A{2}(6,5)); % GPi -> Thal
        X2(6) = E(6).*exp(pinst.A{1}(1,6)); % Thal -> CTX
        
        % X3) compute the "hyperdirect" loop - CTX->STN->GPi->Thal->CTX
        X3(1) = E(1).*exp(pinst.A{1}(4,1)); % M2 -> STN
        X3(2) = E(4).*exp(pinst.A{1}(5,4)); % STN -> GPi
        X3(3) = E(5).*exp(pinst.A{2}(6,5)); % GPi -> Thal
        X3(4) = E(6).*exp(pinst.A{1}(1,6)); % Thal -> CTX
        
        % X4) compute the "thalamocortical" loop - CTX->Thal->CTX
        X4(1) = E(1).*exp(pinst.A{1}(6,1)); % M2 -> STN
        X4(2) = E(6).*exp(pinst.A{1}(1,6)); % STN -> GPi
        
        % X5) compute the "pallido-subthalamic" loop - GPe->STN->GPe
        X5(1) = E(3).*exp(pinst.A{2}(4,3)); % GPe -> STN
        X5(2) = E(4).*exp(pinst.A{1}(3,4)); % STN -> GPe
        
        % Compute Summaries
        netA(:,i) = [sum(X1) sum(X2) sum(X3) sum(X4) sum(X5)]';
        netAbsA(:,i) = [sum(abs(X1)) sum(abs(X2)) sum(abs(X3)) sum(abs(X4)) sum(abs(X5))]';
        
        
        [bpowr_br(i),fpow_br(i),bpowr(i),fpow(i),bcohr(i),fcoh(i)] = computeBetaSpectralStats(R.frqz,{feat});
        if log10(bpowr(i))> -12
            netA(:,i) = nan(1,5);
            netAbsA(:,i) = nan(1,5);
            bpowr(i) = nan;
            bpowr_br(i) = nan;
            bcohr(i) = nan;
            fpow(i) = nan;
        end
    else
        netA(:,i) = nan(1,5);
        netAbsA(:,i) = nan(1,5);
        bpowr(i) = nan;
        bpowr_br(i) = nan;
        bcohr(i) = nan;
        fpow(i) = nan;
    end
end


% Plot BarPlots
% Percentage Diff
ptits = {'% Change in Net Strength','Change in E/I Balance'};
cmap = brewermap(128,'RdBu');
loopNames  = {'Direct','Indirect','Hyperdirect','Thalamocortical','Pallido-subthalamic'};
starheight = [30 1700];
ylax = [-35 35; -2000 2000];
xlax2 = [-50 50; -2500 2500];
for  PL = 1:2
    lowInd = bpowr<prctile(bpowr,25);
    highInd = bpowr>=prctile(bpowr,75);
    
    if PL == 1
        Y = netAbsA;
    else
        Y = netA;
    end
    
    for i = 1:5
        netAXbar = nanmean(Y(i,:));
        if PL == 1
            Z(i,:) = 100.*(Y(i,:)-netAXbar)./netAXbar;
        else
            Z(i,:) = Y(i,:)-netAXbar;
        end
        
        Zlow = Z(i,lowInd);
        Zhigh = Z(i,highInd);
        
        XG(i,1) = mean(Zlow);
        XG(i,2) = mean(Zhigh);
        
        Xvar(i,1) = std(Zlow)./sqrt(numel(Zlow));
        Xvar(i,2) = std(Zhigh)./sqrt(numel(Zhigh));
        
        [h p(i)] = ttest2(Zlow,Zhigh)
    end
    figure(PL)
    B = bar(1:5,XG,'EdgeColor','none');
    B(1).FaceColor = cmap(18,:);
    B(2).FaceColor = cmap(end-18,:);
    hold on
    ax = gca;
    ax.XTickLabel = loopNames; ax.XTickLabelRotation = 45;
    errorbar((1:5)-0.15,XG(:,1),Xvar(:,1),'LineStyle','none','Color','k')
    errorbar((1:5)+0.15,XG(:,2),Xvar(:,2),'LineStyle','none','Color','k')
    
    ylabel(ptits{PL})
    ind = find(p<=0.05 & p>0.01);
    for i = ind
        text(i,starheight(PL),'*','FontSize',18)
    end
    
    ind = find(p<=0.01 & p>0.001);
    for i = ind
        text(i,starheight(PL),'**','FontSize',18)
    end
    
    ind = find(p<=0.001);
    for i = ind
        text(i,starheight(PL),'***','FontSize',18)
    end
    legend(B,{'B Power < 25th','B Power > 25th'},'Location','southeast')
    ylim(ylax(PL,:))
    
    % Now Plot!
    % Plot Spectra
    % for i = 1:100
    %     plotABCSpectraOnly(R.data.feat_xscale,R.data.feat_emp,permMod.feat_rep{i})
    % end
    %
    % for p = 1:4; subplot(1,6,p); ylim([0 1]); end
    
    % Plot Scatters
    figure(PL*10)
    
    for i = 1:5
        subplot(1,5,i)
        scatter((Z(i,:)),log10(bpowr(:)),'MarkerEdgeColor',[0.15 0.15 0.15]);
        hold on
        scatter((Z(i,lowInd)),log10(bpowr(lowInd)),'MarkerEdgeColor',cmap(18,:),'MarkerFaceColor',cmap(18,:));
        scatter((Z(i,highInd)),log10(bpowr(highInd)),'MarkerEdgeColor',cmap(end-18,:),'MarkerFaceColor',cmap(end-18,:))
        
        xlabel(ptits(PL)); ylabel('log Beta Power')
        title(loopNames{i})
        xlim(xlax2(PL,:)); 
        ylim([-14.5 -12])
    end
        set(gcf,'Position',[26         714        1638         264]);

end