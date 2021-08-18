function burstPropertiesCheck(R,fresh)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
close all

% If doing cond compar
if fresh == 1
    a= 0;
    for CON = [1 3]; %[1 3]
        a = a+1;
        load([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF.mat'],'BB')
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
                % Retrieve Burst Envelopes
                rawEnvs = BB.segEnv{K};
                % Compute Threshold Exceedence
                try
                    burstDur(:,chlock,K,a) = meanCI(BB.segDur{K});
                    burstAmp(:,chlock,K,a) = meanCI(BB.segAmp{K});
                    burstInt(:,chlock,K,a) = meanCI(BB.segInterval{K});
                    burstRate(:,chlock,K,a) = meanCI(BB.segRate{K});
                catch
                    disp('This condition didnt yield bursts!')
                end
            end
        end
    end
    
    mkdir([R.rootn '\data\Burst_coinc'])
    save([R.rootn '\data\Burst_coinc\burst_prop.mat'],'burstDur','burstAmp','burstInt','burstRate','R')
else
    load([R.rootn '\data\Burst_coinc\burst_prop.mat'],'burstDur','burstAmp','burstInt','burstRate','R')
end

rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];

    bsel = [1 5 10 15 20 25];
cmap = linspecer(5);
    a= 0;
for CON = [1 3]
    a = a+1;
    load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
    propnames = {'burstDur','burstRate'};
    for prop = 1:2
    figure(1)
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
            yl = [-60 60];
        elseif prop == 2
                        yyaxis right
            si = 'square';
            ls = '--';
            yl = [-20 20];
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
%         ylim(yl)
        grid on
        ax = gca;
%         ax.YTick = 0:25:100;
    end
    set(gcf,'Position',[274 230 1119 532])
    end
end

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
