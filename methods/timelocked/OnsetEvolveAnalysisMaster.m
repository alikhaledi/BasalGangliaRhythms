function OnsetEvolveAnalysisMaster(R)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
BB.threshold_type = 'baseModelThresh';
% R.BB.threshold_type = 'localThresh';

% If doing cond compar
condsel = 1:25;
% close all
a= 0;
for CON = [1 3]
    a = a+1;
    load([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF' BB.threshold_type '.mat'],'BB')
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
        TL = [];

    BB.struccmap = linspecer(4);
    TL.struccmap = BB.struccmap;
    
    condOnsetStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    condPeakStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    condTermStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    for cond = condsel
%         TL.periodT = [-250 250 250]; % beginning end ofset end (to run OLD analysis)
        TL.periodT = [-150 150 600]; % beginning end ofset end
        TL = defineBurstTimeLockEpoch2(BB,TL,cond);
        
        
        if isfield(TL,'onsetT') && size(TL.onsetT{cond},2)>3
            Z = 1.96; % Z value for 95% confidence ilevel
            % Get Onset Stats
            condOnsetStat(:,cond,1) = nanmean(TL.onsetT{cond}');
            condOnsetStat(:,cond,2) = Z.*nanstd(TL.onsetT{cond}')./sqrt(sum(~isnan(TL.onsetT{cond}),2))'; %(prctile(TL.onsetT{cond}',66)-prctile(TL.onsetT{cond}',16))./sqrt(size(TL.onsetT{cond},2));
            
            condTermStat(:,cond,1) = nanmean(TL.onsetOffT{cond}');
            condTermStat(:,cond,2) = Z.*nanstd(TL.onsetOffT{cond}')./sqrt(size(TL.onsetOffT{cond},2))'; %(prctile(TL.onsetOffT{cond}',66)-prctile(TL.onsetOffT{cond}',16))./sqrt(size(TL.onsetOffT{cond},2));
            
            condBurstSuccess(:,cond) = 100.*sum(~isnan(TL.onsetT{cond}),2)./size(TL.onsetT{cond},2)';
            
        else
            condOnsetStat(:,cond,1) = nan;
            condOnsetStat(:,cond,2) = nan;
            
            condTermStat(:,cond,1) = nan;
            condTermStat(:,cond,2) = nan;
            
            condBurstSuccess(:,cond) = nan(4,1);
        end
    end
    pOnset = []; pOffset = [];
    for cond = condsel
        for i = 1:3
            if size(TL.onsetT{cond},2)>3
                [pOnset(i,cond)] = ranksum(TL.onsetT{cond}(4,:),TL.onsetT{cond}(i,:));
                [pOffset(i,cond)] = ranksum(TL.onsetOffT{cond}(4,:),TL.onsetOffT{cond}(i,:));
                
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
    
    subplot(1,2,a)
%     betalist = linspace(10,190,19);
    Klist = BB.condlist(CON,:)*100;
    
    
    for i = 1:4
        pinds = find(~isnan(condOnsetStat(i,:,1)));
        [al aa(i)] = boundedline(Klist(pinds),condOnsetStat(i,pinds,1),[condOnsetStat(i,pinds,2); condOnsetStat(i,pinds,2)]');
        al.Color = BB.struccmap(i,:);
        al.LineWidth = 2;
        al.LineStyle = ':';
        aa(i).FaceColor = BB.struccmap(i,:);
        aa(i).FaceAlpha = 0.4;
        
       
        [cl ca(i)] = boundedline(Klist(pinds),condTermStat(i,pinds,1),[condTermStat(i,pinds,2); condTermStat(i,pinds,2)]');
        cl.Color = BB.struccmap(i,:);
        cl.LineStyle = '--';
        cl.LineWidth = 2;
        ca(i).FaceColor = BB.struccmap(i,:);
        ca(i).FaceAlpha = 0.4;
        
        ylim([-100 450]); grid on
        xlabel('Strength Elliciting %Beta')
        ylabel('Onset Time (ms)')
    end
    set(gca, 'YDir', 'reverse')

    if CON == 1
        xlim([0 500])
    elseif CON == 3
        xlim([0 125])
    end
    X = [20 420];
    offset = 5;
    % Now Plot Sig Bars
    for i = 1:3
        %% ONSET STARS
        bX = pOnsetStar(i,:).*X(1);
        bX(bX==0) = nan;
        bX = bX + (i-1)*offset;
        hold on
        cp = plot(Klist,bX,'Color',BB.struccmap(i,:),'LineWidth',2);
        
        %% OFFSET STARS
        bX = pOffsetStar(i,:).*X(2);
        bX(bX==0) = nan;
        bX = bX + (i-1)*offset;
        hold on
        cp = plot(Klist,bX,'Color',BB.struccmap(i,:),'LineWidth',2);
    end
    
end
set(gcf,'Position',[488.0000   66.6000  599.4000  695.4000])
a = 1;


%%% SCRIPT GRAVE
% NP CI Calculation
% function ci = npCI(X)
% N = numel(X);
% lb = (N/2) - ((1.96*sqrt(N))/2);
% ub = 1+ (N/2) + ((1.96*sqrt(N))/2);
% 
% X = sort(X);
% ub = fix(ub); lb = fix(lb);
% if (ub>0) & (lb>0)
%     ub = X(fix(ub));
%     lb = X(fix(lb));
%     ci = ub-lb;
% else
%     ci = nan;
% end