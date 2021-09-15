function R = SI_CrossCorrelation_corticosubthalamic(R)
close all
rootan = [R.rootn 'data\ConnectionSweep'];

R.CONnames = {'M2 -> STN','GPe -| STN'};
R.condname = {'Fitted','HD-Down','HD-Up','Fitted','PS-Down','PS-Up'};
cmap1 = brewermap(40,'Greys');
cmap2 = brewermap(40,'Blues');
scmap = brewermap(4,'Set1');
statecmap{1} = [0 0 0; scmap(1:2,:)];
statecmap{2} = [0 0 0; scmap(3:4,:)];
hdext = '_REV'; % version tag

for CON = 1:2
    load([rootan '\BB_' R.out.tag '_DiscreteData.mat'],'dataSelect')
    
    for state = 1:3
        X = dataSelect{CON}{state}{1}(4,:);
        Y_M2 = dataSelect{CON}{state}{1}(1,:);
        Y_STN = dataSelect{CON}{state}{1}(3,:);
        
        
        [xcf(:,1),lags] = crosscorr(X,Y_M2,'NumLags',300);
        [xcf(:,2),lags] =crosscorr(X,Y_STN,'NumLags',300);
        lags = lags*R.IntP.dt*1000;
        
        figure(1)
        subplot(2,2,CON)
        plot(lags,xcf(:,1),'Color',statecmap{CON}(state,:),'LineWidth',1.5); hold on
        ylabel('M2/STN XCOR'); xlabel('Lag (ms)'); title(['Modulating ' R.CONnames{CON}]); grid on; box off; ylim([-1 1])
        hold on
        subplot(2,2,CON+2)
        plot(lags,xcf(:,2),'Color',statecmap{CON}(state,:),'LineWidth',1.5); hold on
        ylabel('GPe/STN XCOR'); xlabel('Lag (ms)'); grid on; box off; ylim([-1 1])
        if CON == 1
            legend(R.condname(1:3),'Box','off','Orientation','horizontal')
        else CON == 2
            legend(R.condname(4:6),'Box','off','Orientation','horizontal')
        end
        
    end
end
set(gcf,'Position',[ 680          81        1082         897])
a = 1;