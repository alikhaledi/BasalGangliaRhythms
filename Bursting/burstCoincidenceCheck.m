function burstCoincidenceCheck(R)
rootan = [R.rootn 'data\' R.out.oldtag '\ConnectionSweep'];

% If doing cond compar
condsel = 1:19;
% close all
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
        R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
                BB.segEnv = [];

        BB = defineBetaEvents(R,BB,0);
        for K = 1:numel(BB.segEnv)
            % Retrieve Burst Envelopes
            rawEnvs = BB.segEnv{K};
            % Compute Threshold Exceedence
            try
                binEnv =     rawEnvs>BB.epsAmpfull;
                binEnvSamps = squeeze(sum(binEnv,2));
                binEnvSamps = binEnvSamps(:,3:end);
                % Find bursts that have samples exceeding limit
                cooinc_stat{CON}(:,chlock,K) = sum(binEnvSamps>BB.period,2)./size(binEnvSamps,2);
            catch
                cooinc_stat{CON}(:,chlock,K) = nan(4,1,1);
                disp('This condition didnt yield bursts!')
            end
        end
    end
end


betalist = linspace(10,190,19);
for CON = [1 3]
    
    for chlock = 1:4
        X = squeeze(cooinc_stat{CON}(:,chlock,:));
        if CON == 1
            subplot(2,4,chlock)
        elseif CON == 3
            subplot(2,4,chlock+4)
        end
        plot(betalist,X','LineWidth',2);
        title([R.chloc_name{chlock} '  burst def'])
        if CON == 3
            set(gca, 'XDir', 'reverse')
        end
        xlabel('STN Beta Ellicited %Base')
        ylabel('Percentage of co-incident bursts')
    end
    
end

legend(R.chloc_name)

