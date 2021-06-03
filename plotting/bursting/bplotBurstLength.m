function bplotBurstLength(segA,segL,statename,statecmap,CON)
philist = [100 6 12];

for band = 1:2
    for feat = 1:2
        if feat == 1
            XY = segA;
            ttest2(XY{band,1}{1}(1,:),XY{band,1}{1}(1,:))
            stvec(:,1) = statvec(XY{band,1,CON}{1}(1,:),XY{band,2,CON}{1}(1,:),1);
            stvec(:,2) = statvec(XY{band,1,CON}{1}(1,:),XY{band,3,CON}{1}(1,:),1);
            stvec(:,3) = statvec(XY{band,1,CON}{4}(4,:),XY{band,2,CON}{4}(4,:),1);
            stvec(:,4) = statvec(XY{band,1,CON}{4}(4,:),XY{band,3,CON}{4}(4,:),1);
            N1 = median(XY{band,1,CON}{1}(1,:)); N2 = median(XY{band,1,CON}{4}(4,:));
            
                    X = {perCh(XY{band,2,CON}{1}(1,:),N1); perCh(XY{band,1,CON}{1}(1,:),N1); perCh(XY{band,3,CON}{1}(1,:),N1);...
            perCh(XY{band,1,CON}{4}(4,:),N2); perCh(XY{band,1,CON}{4}(4,:),N2); perCh(XY{band,3,CON}{4}(4,:),N2)};
        elseif feat ==2
            XY = segL;
            ttest2(XY{band,1}{1}(1,:),XY{band,1}{1}(1,:))
            stvec(:,1) = statvec(XY{band,1,CON}{1},XY{band,2,CON}{1},1);
            stvec(:,2) = statvec(XY{band,1,CON}{1},XY{band,3,CON}{1},1);
            stvec(:,3) = statvec(XY{band,1,CON}{4},XY{band,2,CON}{4},1);
            stvec(:,4) = statvec(XY{band,1,CON}{4},XY{band,3,CON}{4},1);
            
           N1 = median(XY{band,1,CON}{1}); N2 = median(XY{band,1,CON}{4});
                    X = {XY{band,2,CON}{1}- N1; XY{band,1,CON}{1}- N1; XY{band,3,CON}{1}- N1;...
            XY{band,2,CON}{4}- N2; XY{band,1,CON}{4}- N2; XY{band,3,CON}{4}- N2};

        end
    subplot(2,2,sub2ind([2 2],feat,band))
        
        G = []; Gc = []; Gcmap = [];
        gc = [2 1 3 2 1 3];
        for g = 1:size(X,1)
            G = [G repmat(g,1,numel(X{g}))];
            Gc = [Gc; repmat(statecmap(gc(g),:),numel(X{g}),1)];
            Gcmap = [Gcmap repmat(gc(g),1,numel(X{g}))]
        end
        
        X = [X{:}];
        boxplot(X,G,'colorgroup',Gcmap,'colors',statecmap,'BoxStyle','filled','Widths',0.5,'symbol','+','jitter',0.5)
        
        a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',20); % Set width
idx=strcmpi(t,'Whisker');  % Find Box objects
set(a(idx),'linewidth',2); % Set width
idx=strcmpi(t,'Median');  % Find Box objects
set(a(idx),'linewidth',2); % Set width
        grid on; box off;
        
        if feat == 1
            ylabel('%Change in burst amplitude')
        elseif feat == 2
            ylabel('Change in burst duration (ms)')
        end
        
        a = gca;
        a.XTickLabel = {[statename '-Down M2'],'Base M2',[statename '-Up M2'],...
            [statename '-Down STN'],'Base STN',[statename '-Up STN']}
        a.XTickLabelRotation = -45;
    end
end

function PC = perCh(X,N)
PC = (X-N)./N;
PC = PC*100;


