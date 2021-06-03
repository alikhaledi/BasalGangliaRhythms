function plotOverlapBars(sampBurstOverl)

scmap = brewermap(4,'Set1');
% scmap = [0 0 0; scmap];
labnames = {'M2/STN','STR/STN','GPe/STN','STN/STN','GPi/STN','Thal./STN'};
titname = {'Change in burst overlap for low beta','Change in burst overlap for high beta'};

    
for CON = [1 3]
    for state = 1:3
        for ch = 1:6
            for band = 1:2
                %             XS(ch,state,CON) = sum(sampBurstOverl{band,state,CON}{ch,4}>0)./numel(sampBurstOverl{band,state,CON}{ch,4});
                XS(band,ch,state,CON) = (nanmean(sampBurstOverl{band,state,CON}{ch,4})-nanmean(sampBurstOverl{band,1,CON}{ch,4}))./nanmean(sampBurstOverl{band,1,CON}{ch,4});
                XH(band,ch,state,CON) = (nanstd(sampBurstOverl{band,state,CON}{ch,4})-nanstd(sampBurstOverl{band,1,CON}{ch,4}))./nanstd(sampBurstOverl{band,1,CON}{ch,4});
            end
        end
    end
end


    YS = [];
    YS = [YS; squeeze(XS(1,:,2,1)); squeeze(XS(2,:,2,1)); squeeze(XS(1,:,3,1)); squeeze(XS(2,:,3,1))]; % CON = 1
    YS = [YS; squeeze(XS(1,:,2,3)); squeeze(XS(2,:,2,3)); squeeze(XS(1,:,3,3)); squeeze(XS(2,:,3,3))]; % CON = 3
    
    %     p(:,:,2) = [];
%     p = [p(:,2:3,1) p(:,2:3,2)];
    
    b = bar(YS');
    for ib =1:4
        b1 = ((2*ib)-1);
        b2 = (2*ib);
        b(b1).FaceColor = scmap(ib,:);
b(b2).FaceColor = scmap(ib,:).*0.7;
    end
    
    a = gca;
    a.XTickLabel =     labnames;
    a.XTickLabelRotation = -45;
    grid on; box off
    title(titname{band})
    ylabel('% Change from baseline')
legend({'HD-Down - LB','HD-Down - HB','HD-Up - LB','HD-Up - HB','PS-Down - LB','PS-Down - HB','PS-Up - LB','PS-Up - HB'})
