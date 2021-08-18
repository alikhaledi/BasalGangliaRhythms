function burstCoincidenceCheck(R,fresh)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
close all
% BB.threshold_type = 'localThresh';
R.BB.threshold_type = 'baseModelThresh';

% If doing cond compar
if fresh == 1
    a= 0;
    for CON = [1 3]; %[1 3]
        a = a+1;
        load([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF' R.BB.threshold_type '.mat'],'BB')
        load([rootan '\BB_'  R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
        BB.struccmap = linspecer(4);
        TL.struccmap = BB.struccmap;
        for chlock = 1:4 % loop through burst channel def;
            % Recompute Bursts changing the channel lock definition
            R1 = R;
            R1.condname = {};
            R = compute_BetaBursts_Simulated(R1);
            R.BB.pairInd = [1 chlock];
            R.BB.minBBlength = 1; %  Minimum burst period- cycles
            BB.segEnv = [];
            
            BB = defineBetaEvents(R,BB,0);
            for K = 1:numel(BB.segEnv)
                if numel(BB.segEnv{K})>0
                    % Retrieve Burst Envelopes
                    rawEnvs = BB.segEnv{K};
%                     if chlock == 4
%                         try
%                             Cvar = nanstd(reshape(BB.segEnv{K},4,[]),[],2);
%                             BB.segEnv{K} = (BB.segEnv{K})./Cvar;
%                             rawEnvs(:,:,1,K,a)= squeeze(mean(BB.segEnv{K},3));
%                             rawEnvs(:,:,2,K,a)= Z.*squeeze(std(BB.segEnv{K},[],3))./sqrt(size(BB.segEnv{K},3));
%                         catch
%                             rawEnvs(:,:,1,K,a)= nan(4,size(rawEnvs,2));
%                             rawEnvs(:,:,2,K,a)= nan(4,size(rawEnvs,2));
%                         end
%                     end
                    % Compute Threshold Exceedence
                    switch BB.threshold_type
                        case 'localThresh'
                            binEnv =     rawEnvs>BB.epsAmpfull(:,K);
                        case 'baseModelThresh'
                            binEnv =     rawEnvs>BB.epsAmpfull(:,1);
                    end
                    
                    burstBool = [];
                    for seg = 3:size(binEnv,3)
                        for ch = 1:4
                            X = find(squeeze(binEnv(ch,:,seg)));
                            L = SplitVec(X,'consecutive');
                            segL = cellfun('length',L);
                            ab = find(segL>BB.period);
                            
                            %                             if numel(ab)>1
                            %                                 for i = 1:numel(ab)
                            %                                     seginit(i) = L{i}(1);
                            %                                 end
                            %                                 [dum segch] = min(abs(seginit-(size(tvec,2)/2)));
                            %                             else
                            %                                 segch = 1;
                            %                             end
                            
                            if ~isempty(ab)
                                burstBool(ch,seg)= 1;
                            else
                                burstBool(ch,seg)= 0;
                            end
                        end                       
                    end
                    burstBool(:,1:2) = [];
                    
                    % Find bursts that have samples exceeding limit
                    cooinc_stat{CON}(:,chlock,K) = sum(burstBool,2)./size(burstBool,2);
                    burstDur(:,chlock,K,a) = meanCI(BB.segDur{K});
                    burstAmp(:,chlock,K,a) = meanCI(BB.segAmp{K});
                    burstInt(:,chlock,K,a) = meanCI(BB.segInterval{K});
                    burstRate(:,chlock,K,a) = meanCI(BB.segRate{K});
                else
                    cooinc_stat{CON}(:,chlock,K) = nan;
                    burstDur(:,chlock,K,a) = nan;
                    burstAmp(:,chlock,K,a) = nan;
                    burstInt(:,chlock,K,a) = nan;
                    burstRate(:,chlock,K,a) = nan;
                end
            end
        end
    end
    

            TWin = linspace(-150,150,size(rawEnvs,2));

    mkdir([R.rootn '\data\Burst_coinc'])
    save([R.rootn '\data\Burst_coinc\burst_coinc' R.BB.threshold_type '.mat'],'burstDur','burstAmp','burstInt','burstDur','burstRate','TWin','rawEnvs','cooinc_stat','R')
else
    load([R.rootn '\data\Burst_coinc\burst_coinc' R.BB.threshold_type '.mat'],'burstDur','burstAmp','burstInt','burstDur','burstRate','TWin','rawEnvs','cooinc_stat','R')
end


%% PLOT COINCIDENCE
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
    bsel = [1 5 10 15 21 25];
figure
cmap = linspecer(5);
for CON = [1 3]
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
    for i = 1:4
        if CON == 1 
            subplot(2,4,i)
            %             K = repmat([10 50 100 150 190],4,1);
            Z = repmat(ck_1(1,:)*100,4,1);
            K = repmat((ck_1(1,bsel))*100,4,1);
        else
            subplot(2,4,i+4)
            %             K = repmat([10 50 100 150 190],4,1);
            Z = repmat(ck_1(3,:)*100,4,1);
            K = repmat((ck_1(3,bsel))*100,4,1);
        end
        X = squeeze(cooinc_stat{CON}(:,i,bsel))*100;
        cind = repmat(1:4,numel(bsel),1)';
        scatter(K(:),X(:),70,cmap(cind(:),:),'filled'); %X(:).*0.7
        hold on
        %         X = squeeze(cooinc_stat{CON}(:,i,:))*100;
        p = plot(K',X','LineWidth',1.5); %X(:).*0.7
        for l = 1:4
            p(l).Color = cmap(l,:);
        end
        xlabel('Connection %')
        ylabel('% Burst Coincidence')
        title(['Locked to ' R.chloc_name{i} ' bursts'])
        if CON == 1
            xlim([0 500])
        elseif CON == 3
            xlim([0 125])
        end
        ylim([0 125])
        grid on
        a = gca;
        a.YTick = 0:25:100;
    end
    
    set(gcf,'Position',[274 230 1119 532])
end

%% PLOT PROPERTIES
figure

    bsel = [1 5 10 15 21 25];
cmap = linspecer(5);
    a= 0;
for CON = [1 3]
    a = a+1;
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
    propnames = {'burstDur','burstRate'};
    for prop = 1:2
    for i = 1:4
        if CON == 1
            subplot(2,4,i)
%             K = repmat([10 50 100 150 190],1,1);
            K = repmat((ck_1(1,bsel))*100,4,1);
        else
            subplot(2,4,i+4)
%             K = repmat([10 50 100 150 190],1,1);
            K = repmat((ck_1(3,bsel))*100,4,1);
        end
        
        if prop == 1
            yyaxis left
            si = 'o';
            ls = '-';
            yl = [0 1200];
        elseif prop == 2
                        yyaxis right
            si = 'square';
            ls = '--';
            yl = [0 200];
        end
        
        burstprop = eval(propnames{prop});
        X = squeeze(burstprop(1,i,bsel,a));
%         X = ((X-X(3))./X(3)).*100;
        
        Xv = squeeze(burstprop(2,i,bsel,a));
%         Xv = ((Xv-Xv(3))./Xv(3)).*100;
        
        scatter(K(1,1:6),X,70,cmap(i,:),'filled','Marker',si); %X(:).*0.7
        hold on
        p = plot(K(1,1:6),X); %X(:).*0.7
        p.Color = cmap(i,:);
        p.LineStyle = ls;

        xlabel('% Connection Strength')
        ylabel([propnames{prop} '(ms)'])
        title(['Bursts defined at ' R.chloc_name{i}])
%         xlim([0 200])
        ylim(yl)
        grid on
        ax = gca;
                if CON == 1
            xlim([0 500])
        elseif CON == 3
            xlim([0 125])
        end
%         ax.YTick = 0:25:100;
    end
    set(gcf,'Position',[274 230 1119 532])
    end
end

%% PLOT LOCKANALY
figure

bEquiv = {'10%','100%','Max%'};
for con = 1:2
    if con == 1
        CON = 1;
    else
        CON = 3;
    end
    bsel = [];
            bsel{CON} = [1 find(ck_1(CON,:)==1) 25];

    p = 0;
    for K = bsel{con}
        p = p+1;
        subplot(2,3,sub2ind([3 2],p,con))
        for ch = 1:4
            [ha(ch) hb] = boundedline(TWin,squeeze(rawEnvs(ch,:,1,K,con)),squeeze(rawEnvs(ch,:,2,K,con)))
            ha(ch).Color = cmap(ch,:);
            ha(ch).LineWidth = 2;
            hb.FaceColor = cmap(ch,:);
            hb.FaceAlpha = 0.7;
            hold on
%             if ch ==4
%                 plot([-150 150],[ epsAmp{con}(ch)  epsAmp{con}(ch)],'Color',cmap(ch,:),'LineStyle','--','LineWidth',2)
%             end
        end
        xlabel('Time relative to STN Onset')
        ylabel('Beta Envelope')
        title(['Connection ' bEquiv{p}])
        %             xlim([0 200])
        %             ylim([0 125])
        grid on
        a = gca;
    end
    legend(ha,R.chloc_name)
end
set(gcf,'Position',[274 230 1119 532])



function stat = meanCI(X)
 stat(1) = nanmean(X);
Z = 1.96;
 stat(2) = Z.*squeeze(nanstd(X))./sqrt(sum(~isnan(X)));
%% SCRIPT GRAVE
% Old Plotting
% for CON = [1 3]
%
%     for chlock = 1:4
%         X = squeeze(cooinc_stat{CON}(:,chlock,:));
%         if CON == 1
%             subplot(2,4,chlock)
%         elseif CON == 3
%             subplot(2,4,chlock+4)
%         end
%         plot(betalist,X','LineWidth',2);
%         title([R.chloc_name{chlock} '  burst def'])
%         if CON == 3
%             set(gca, 'XDir', 'reverse')
%         end
%         xlabel('STN Beta Ellicited %Base')
%         ylabel('Percentage of co-incident bursts')
%     end
%
% end
%
% legend(R.chloc_name)
%
