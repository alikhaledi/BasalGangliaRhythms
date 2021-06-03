function plotOverlapBars(overlapMat)
scmap = brewermap(4,'Set1');
labnames = {'M2/STN','STR/STN','GPe/STN','GPi/STN','Thal./STN'};
titname = {'Change in burst overlap for low beta','Change in burst overlap for high beta'};
for band = 1:2
    
    for CON = [1 3]
        for state = 1:3
            
            X1 = [overlapMat{band,state,CON}(:,4,1); overlapMat{band,state,CON}(4,:,1)'];
            X1(isnan(X1) | X1==0) = [];
            
            S1 = [overlapMat{band,state,CON}(:,4,2); overlapMat{band,state,CON}(4,:,2)'];
            S1(isnan(S1) | S1==0) = [];
            
            N1 = [overlapMat{band,state,CON}(:,4,3); overlapMat{band,state,CON}(4,:,3)'];
            N1(isnan(N1) | N1==0) = [];
            
            
            X2 = [overlapMat{band,1,CON}(:,4,1); overlapMat{band,1,CON}(4,:,1)'];
            X2(isnan(X2) | X2==0) = [];
            
            S2 = [overlapMat{band,1,CON}(:,4,2); overlapMat{band,1,CON}(4,:,2)'];
            S2(isnan(S2) | S2==0) = [];
            
            N2 = [overlapMat{band,1,CON}(:,4,3); overlapMat{band,1,CON}(4,:,3)'];
            N2(isnan(N2) | N2==0) = [];
            
            p(:,state,CON)  = indepTTest(X1,X2,S1,S2,N1,N2);
            
            XS(:,state,CON) = X1;
            
        end
    end
    
    XS(:,:,2) = [];
    XS = [XS(:,1:3,1) XS(:,1:3,2)];
    
    p(:,:,2) = [];
    p = [p(:,2:3,1) p(:,2:3,2)];
    
    subplot(2,1,band)
    b = bar(XS);
    for ib =1:4
        b(ib).FaceColor = scmap(ib,:);
    end
    
    a = gca;
    a.XTickLabel =     labnames;
    a.XTickLabelRotation = -45;
    grid on; box off
    title(titname{band})
    ylabel('% Change from baseline')
end
