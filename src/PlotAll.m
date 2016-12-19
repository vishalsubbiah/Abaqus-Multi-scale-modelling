load('FinResults1-10')
load('FinalResults1-10')

h10=figure();
scatter(Results(:,1)*180/pi,Results(:,2));
set(gca,'FontSize',16);
ylabel('Energy Gradient (J)');
xlabel(' Angle');
saveas(h10,'ScatterPlot10.png');
close(h10)

g10=figure();
plot(FinResults(:,1)*180/pi,FinResults(:,2));
set(gca,'FontSize',16);
ylabel('Energy Gradient (J)');
xlabel('Misorientation Angle');
saveas(g10,'Plot10.png');
close(g10)


load('FinResults1')
load('FinalResults1')

h=figure();
scatter(Results(:,1)*180/pi,Results(:,2));
set(gca,'FontSize',16);
ylabel('Energy Gradient (J)');
xlabel(' Angle');
saveas(h,'ScatterPlot.png');
close(h)

g=figure();
plot(FinResults(:,1)*180/pi,FinResults(:,2));
set(gca,'FontSize',16);
ylabel('Energy Gradient (J)');
xlabel('Misorientation Angle');
saveas(g,'Plot.png');
close(g)

