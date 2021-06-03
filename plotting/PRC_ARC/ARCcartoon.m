cmap = brewermap(12,'Set1');
phangle = [6 12];

phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12


xt = linspace(0,2*pi,130);


plot(rad2deg(xt),sin(xt),'k-','LineWidth',2);
hold on
plot(rad2deg(xt(51:60)),sin(xt(51:60)),'k-','LineWidth',4,'Color',cmap(6-1,:));
hold on
plot(rad2deg(xt(116:125)),sin(xt(116:125)),'k-','LineWidth',4,'Color',cmap(12-1,:));

a = gca;
a.XTick = [0 180 360];
ylim([-1.5 1.5])