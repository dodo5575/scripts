function CurrentAnalysis(name)

TitleSize=16;
LabelSize=14;
AxisSize=12;


%load input data
currentFile1=sprintf('curr_%s.dat', name);
File1=load(currentFile1);

currentFile2=sprintf('currCl_%s.dat', name);
File2=load(currentFile2);

currentFile3=sprintf('currMg_%s.dat', name);
File3=load(currentFile3);

currentFile4=sprintf('currK_%s.dat', name);
File4=load(currentFile4);


currentTraceOut=sprintf('%s_currentTrace', name);
Charge_vs_TimeOut=sprintf('%s_Charge_vs_Time', name);


%extract variables
Time = File1(:,1);

Total = File1(:,2);

Cl = File2(:,2);

Mg = File3(:,2);

K = File4(:,2);


%find the common x and y range for plotting
MinTime = min(Time);
MaxTime = max(Time);
MinX = floor(MinTime/10) * 10;
MaxX = ceil(MaxTime/10) * 10;
Min = min([min(Total),min(Cl),min(Mg),min(K)]);
Max = max([max(Total),max(Cl),max(Mg),max(K)]);
MinY = floor(Min/20) * 20;
MaxY = ceil(Max/20) * 20;



%generate plot for current trace
Figure1=figure;
subplot(2,2,1);
plot(Time,Total,'Color','blue');
title('Total','FontSize',TitleSize);
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Current(nA)','FontSize',LabelSize);
axis([MinX MaxX MinY MaxY]);
set(gca,'FontSize',AxisSize);
set(gca,'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'YTick',MinY:20:MaxY);
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 1 0]);
%et(xlabh,'Position',get(xlabh,'Position') + [0 10 0]);

subplot(2,2,2);
plot(Time,Cl,'Color','blue');
title('Cl','FontSize',TitleSize);
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Current(nA)','FontSize',LabelSize);
axis([MinX MaxX MinY MaxY]);
set(gca,'FontSize',AxisSize);
set(gca,'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'YTick',MinY:20:MaxY);
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 1 0]);
%set(xlabh,'Position',get(xlabh,'Position') + [0 10 0]);

subplot(2,2,3);
plot(Time,Mg,'Color','blue');
title('Mg','FontSize',TitleSize);
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Current(nA)','FontSize',LabelSize);
axis([MinX MaxX MinY MaxY]);
set(gca,'FontSize',AxisSize);
set(gca,'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'YTick',MinY:20:MaxY);
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 1 0]);
%set(xlabh,'Position',get(xlabh,'Position') + [0 10 0]);

subplot(2,2,4);
plot(Time,K,'Color','blue');
title('K','FontSize',TitleSize);
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Current(nA)','FontSize',LabelSize);
axis([MinX MaxX MinY MaxY]);
set(gca,'FontSize',AxisSize);
set(gca,'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'YTick',MinY:20:MaxY);
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 1 0]);
%set(xlabh,'Position',get(xlabh,'Position') + [0 10 0]);
saveas(Figure1,currentTraceOut,'epsc2');


%calculate charge
cumTotal = cumsum(Total)/1000;  
cumCl = cumsum(Cl)/1000;
cumMg = cumsum(Mg)/1000;
cumK = cumsum(K)/1000;


%plot charge vs time
Figure2=figure;
plot(Time,cumTotal,'LineWidth',1);
hold all;
plot(Time,cumCl,'LineWidth',1);
plot(Time,cumMg,'LineWidth',1);
plot(Time,cumK,'LineWidth',1);
hold off;
title('Charge vs Time','FontSize',TitleSize);
legend('Total','Cl','Mg','K','Location','NorthWest');
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Charge(femto Coulomb)','FontSize',LabelSize);
set(gca,'FontSize',AxisSize);
%xlabh = get(gca,'XLabel');
%set(xlabh,'Position',get(xlabh,'Position') - [0 0.1 0]);
%set(xlabh,'Position',get(xlabh,'Position') + [0 10 0]);
saveas(Figure2,Charge_vs_TimeOut,'epsc2')








