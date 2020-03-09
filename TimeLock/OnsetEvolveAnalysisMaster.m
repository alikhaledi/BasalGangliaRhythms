function OnsetEvolveAnalysisMaster(R)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];

% If doing cond compar
condsel = 1:19;
% close all
a= 0;
for CON = [3 1]
    a = a+1;
    load([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF.mat'],'BB')
    BB.struccmap = linspecer(4);
    TL.struccmap = BB.struccmap;
    
    condOnsetStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    condPeakStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    condTermStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    for cond = condsel
        TL.periodT = [-250 250];
        %     TL.periodT = [-50 300];
        TL = defineBurstTimeLockEpoch(BB,TL,cond);
        
        
        if size(TL.onsetT{cond},2)>3
            % Get Onset Stats
            condOnsetStat(:,cond,1) = nanmedian(TL.onsetT{cond}');
            condOnsetStat(:,cond,2) = (prctile(TL.onsetT{cond}',66)-prctile(TL.onsetT{cond}',16))./sqrt(size(TL.onsetT{cond},2));
            
            condPeakStat(:,cond,1) = nanmedian(TL.maxT{cond}');
            condPeakStat(:,cond,2) = (prctile(TL.maxT{cond}',66)-prctile(TL.maxT{cond}',16))./sqrt(size(TL.maxT{cond},2));
            
            condTermStat(:,cond,1) = nanmedian(TL.onsetOffT{cond}');
            condTermStat(:,cond,2) = (prctile(TL.onsetOffT{cond}',66)-prctile(TL.onsetOffT{cond}',16))./sqrt(size(TL.onsetOffT{cond},2));
            
            condBurstSuccess(:,cond) = 100.*sum(~isnan(TL.onsetT{cond}),2)./size(TL.onsetT{cond},2);
            
        else
            condOnsetStat(:,cond,1) = nan;
            condOnsetStat(:,cond,2) = nan;
            
            condPeakStat(:,cond,1) = nan;
            condPeakStat(:,cond,2) = nan;
            
            condTermStat(:,cond,1) = nan;
            condTermStat(:,cond,2) = nan;
            
            condBurstSuccess(:,cond) = nan(4,1);
        end
    end
    
    for cond = condsel
        for i = 1:3
            if size(TL.onsetT{cond},2)>3
                [h pOnset(i,cond)] = ttest(TL.onsetT{cond}(4,:),TL.onsetT{cond}(i,:));
                [h pOffset(i,cond)] = ttest(TL.onsetOffT{cond}(4,:),TL.onsetOffT{cond}(i,:));
                
                
                % Segment Duration
                [Rk pdum] = corrcoef(TL.segDur{cond}(4,:),TL.segDur{cond}(i,:),'rows','complete');
                 pCorrDurr(i,cond) = pdum(2);
                 
                [Rk pdum] = corrcoef(TL.onsetT{cond}(4,:),TL.onsetT{cond}(i,:),'rows','complete');
                pCorrOnset(i,cond) = pdum(2);
                [Rk pdum] = corrcoef(TL.onsetOffT{cond}(4,:),TL.onsetOffT{cond}(i,:),'rows','complete');
                pCorrOffset(i,cond) = pdum(2);
            else
                pOnset(i,cond) = nan;
                pOffset(i,cond) = nan;
                pCorrOnset(i,cond) = nan;
                pCorrOffset(i,cond) = nan;
            end
        end
    end
    
    pOnsetStar = benHochFWER(pOnset,0.01);
    pOffsetStar = benHochFWER(pOffset,0.01);
    
    pCorrOnsetStar =  benHochFWER(pCorrOnset,0.05);
    pCorrOffsetStar =  benHochFWER(pCorrOffset,0.05);
    
    %     condOnsetStat(isnan(condOnsetStat(:,:,:))) = [];
    subplot(1,2,a)
    betalist = linspace(10,190,19);
    for i = 1:4
        pinds = find(~isnan(condOnsetStat(i,:,1)));
        [al aa(i)] = boundedline(betalist(pinds),condOnsetStat(i,pinds,1),[condOnsetStat(i,pinds,2); condOnsetStat(i,pinds,2)]');
        al.Color = BB.struccmap(i,:);
        al.LineWidth = 2;
        al.LineStyle = ':';
        aa(i).FaceColor = BB.struccmap(i,:);
        aa(i).FaceAlpha = 0.6;
        
        %         [bl ba(i)] = boundedline(betalist,condPeakStat(i,:,1),[condPeakStat(i,:,2); condPeakStat(i,:,2)]')
        %         bl.Color = BB.struccmap(i,:);
        %         bl.LineWidth = 2;
        %         bl.LineStyle = '-';
        %         ba(i).FaceColor = BB.struccmap(i,:);
        %         ba(i).FaceAlpha = 0.6;
        
        [cl ca(i)] = boundedline(betalist(pinds),condTermStat(i,pinds,1),[condTermStat(i,pinds,2); condTermStat(i,pinds,2)]');
        cl.Color = BB.struccmap(i,:);
        cl.LineStyle = '--';
        cl.LineWidth = 2;
        ca(i).FaceColor = BB.struccmap(i,:);
        ca(i).FaceAlpha = 0.6;
        
        %         plot(betalist,condPeakStat(i,:,1),'color',BB.struccmap(i,:),'LineWidth',2)
        %         hold on
        %         plot(betalist,condOnsetStat(i,:,1),'color',BB.struccmap(i,:),'LineWidth',2,'LineStyle',':')
        %         plot(betalist,condTermStat(i,:,1),'color',BB.struccmap(i,:),'LineWidth',2,'LineStyle','--')
        
        xlim([10 190]); ylim([-75 130]); grid on
        xlabel('Strength Elliciting %Beta')
        ylabel('Onset Time (ms)')
    end
    set(gca, 'YDir', 'reverse')
    if CON == 1
        set(gca, 'XDir', 'reverse')
    end
    X = [5 110];
    offset = 2.5;
    % Now Plot Sig Bars
    for i = 1:3
        %% ONSET STARS
        bX = pOnsetStar(i,:).*X(1);
        bX(bX==0) = nan;
        bX = bX + (i-1)*offset;
        hold on
        cp = plot(betalist,bX,'Color',BB.struccmap(i,:),'LineWidth',2);
        %         cp.MarkerFaceColor = BB.struccmap(i,:);
        %         cp.Marker = 'o';
        
% %         bX = pCorrOnsetStar(i,:).*X(1);
% %         bX(bX==0) = nan;
% %         bX = bX + (i-1)*offset;
% %         hold on
% %         csp = scatter(betalist,bX);
% %         csp.Marker = 'o';
% %         csp.MarkerFaceColor = BB.struccmap(i,:);
% %         csp.MarkerEdgeColor = 'none';
% %         csp.SizeData = 25;
        
        
        
        %% OFFSET STARS
        bX = pOffsetStar(i,:).*X(2);
        bX(bX==0) = nan;
        bX = bX + (i-1)*offset;
        hold on
        cp = plot(betalist,bX,'Color',BB.struccmap(i,:),'LineWidth',2);
        %         cp.MarkerFaceColor = BB.struccmap(i,:);
        %         cp.Marker = 'o';
        
        bX = pCorrOffsetStar(i,:).*X(2);
        bX(bX==0) = nan;
        bX = bX + (i-1)*offset;
        hold on
        csp = scatter(betalist,bX);
        csp.Marker = 'o';
        csp.MarkerFaceColor = BB.struccmap(i,:);
        csp.MarkerEdgeColor = 'none';
        csp.SizeData = 25;
    end
    
end
