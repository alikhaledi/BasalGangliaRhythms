function burstLockCheck(R,fresh)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];
R.BB.threshold_type = 'baseModelThresh';

close all
cmap = linspecer(4);
if fresh == 1
    % If doing cond compar
    a= 0;
    for CON = [1 3]; %[1 3]
        % These are the connection weight definitions
        a = a+1;
        load([rootan '\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF.mat'],'BB')
        load([rootan '\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1_bKF.mat'],'ck_1')
        bsel{a} = [1 find(ck_1(CON,:)==1) 25];
        
        
        
        BB.struccmap = linspecer(4);
        TL.struccmap = BB.struccmap;
        % Recompute Bursts changing the channel lock definition
        R1 = R;
        R1.condname = {};
        R = compute_BetaBursts_Simulated(R1);
        R.BB.pairInd = [1 4];
        R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
        BB.segEnv = [];
        
        BB = defineBetaEvents(R,BB,0);
        Z = 1.96; %95% confidence interval
        
        for K = bsel{a}
            % Retrieve Burst Envelopes
            try
                rawEnvs_tmp = BB.segEnv{K};
                switch BB.threshold_type
                    case 'localThresh'
                        binEnv =     rawEnvs_tmp>BB.epsAmpfull(:,K);
                    case 'baseModelThresh'
                        binEnv =     rawEnvs_tmp>BB.epsAmpfull(:,1);
                end
                
                burstBool = [];
                for sch = 1:2
                    for seg = 3:size(binEnv,3)
                        X = find(squeeze(binEnv(sch,:,seg)));
                        L = SplitVec(X,'consecutive');
                        segL = cellfun('length',L);
                        ab = find(segL>BB.period);
                        
                        if ~isempty(ab)
                            burstBool(sch,seg)= 1;
                        else
                            burstBool(sch,seg)= 0;
                        end
                    end
                end
                segstore = []; as = 0;
                for seg = 3:size(BB.segEnv{K},3)
                    if burstBool(2,seg) %and(~burstBool(1,seg),~burstBool(2,seg))
                        as = as+1;
                        X = squeeze(BB.segEnv{K}(:,:,seg));
                        X = (X-mean(X,2))./std(X,[],2);
                        segstore(:,:,as) = X;
                    end
                end
                %                 Cvar = nanstd(reshape(BB.segEnv{K},4,[]),[],2);
                %                 BB.segEnv{K} = (BB.segEnv{K})./Cvar;
                rawEnvs(:,:,1,K,a)= squeeze(mean(segstore,3));
                rawEnvs(:,:,2,K,a)= Z.*squeeze(std(segstore,[],3))./sqrt(size(segstore,3));
            catch
                rawEnvs(:,:,1,K,a)= nan(4,size(rawEnvs,2));
                rawEnvs(:,:,2,K,a)= nan(4,size(rawEnvs,2));
            end
        end
        TWin = linspace(-150,150,size(rawEnvs,2));
        epsAmp{a} = BB.epsAmpfull; %./Cvar;
    end
    mkdir([R.rootn '\data\Burst_coinc'])
    save([R.rootn '\data\Burst_coinc\burstLockEnvs.mat'],'TWin','rawEnvs','epsAmp','bsel')
else
    mkdir([R.rootn '\data\Burst_coinc'])
    load([R.rootn '\data\Burst_coinc\burstLockEnvs.mat'],'TWin','rawEnvs','epsAmp','bsel')
end

%% Plot all 4 Per Conditions
bEquiv = {'10%','100%','Max%'};
for con = 1:2
    
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


%% Plot all 4 Per Conditions
% bsel = [1 10 19];
% bEquiv = {'10%','100%','190%'};
% for con = 1:2
%     p = 0;
%     for K = bsel
%         p = p+1;
%         subplot(2,3,sub2ind([3 2],p,con))
%         for ch = 1:4
%             [ha(ch) hb] = boundedline(TWin,squeeze(rawEnvs(ch,:,1,K,con)),squeeze(rawEnvs(ch,:,2,K,con)))
%             ha(ch).Color = cmap(ch,:);
%             ha(ch).LineWidth = 2;
%             hb.FaceColor = cmap(ch,:);
%             hb.FaceAlpha = 0.7;
%             hold on
%             plot([-150 150],[ epsAmp{con}(ch)  epsAmp{con}(ch)],'Color',cmap(ch,:),'LineStyle','--','LineWidth',2)
%         end
%         xlabel('Time relative to STN Onset')
%         ylabel('% Burst Coincidence')
%         title(['Equivalent Beta power: ' bEquiv{p}])
%         %             xlim([0 200])
%         %             ylim([0 125])
%         grid on
%         a = gca;
%     end
%     legend(ha,R.chloc_name)
% end
% set(gcf,'Position',[274 230 1119 532])



