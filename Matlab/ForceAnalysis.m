function ForceAnalysis(name)

TitleSize=16;
LabelSize=14;
AxisSize=12;


%load input data
forceFile=sprintf('%s_force.dat', name);
File=load(forceFile);



%set output name
forceTraceOut=sprintf('%s_forceTrace', name);
integrativeForceOut=sprintf('%s_momentum', name);


%extract variables
Time = File(:,1);
Force = File(:,2);


%generate plot for force trace
Figure1 = figure;
plot(Time,Force,'Color','Blue');
title('Force Trace','FontSize',TitleSize);
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Force(pN)','FontSize',LabelSize);
saveas(Figure1,forceTraceOut,'epsc2');


%calculate the integrative force
cumForce = cumsum(Force);


%plot integrative force
Figure2 = figure;
plot(Time,cumForce,'Color','Blue','LineWidth',1);
title('Momentum vs Time','FontSize',TitleSize);
xlabel('Time(ns)','FontSize',LabelSize);
ylabel('Momentum(pN ns)','FontSize',LabelSize);
saveas(Figure2,integrativeForceOut,'epsc2')
