function IonProfileAnalysis(name)

TitleSize=32;
LabelSize=28;
AxisSize=24;


%load input data
currentFile1=sprintf('%s_P_profile.dat', name);
File1=load(currentFile1);

currentFile2=sprintf('%s_Cl_profile.dat', name);
File2=load(currentFile2);

currentFile3=sprintf('%s_Mg_profile.dat', name);
File3=load(currentFile3);

currentFile4=sprintf('%s_K_profile.dat', name);
File4=load(currentFile4);


ionProfile=sprintf('%s_IonProfile', name);

%extract variables
zAxis = File1(:,1);

P = File1(:,2);

Cl = File2(:,2);

Mg = File3(:,2);

K = File4(:,2);


%find the common x and y range for plotting
MinZAxis = min(zAxis);
MaxZAxis = max(zAxis);
MinX = floor(MinZAxis/20) * 20;
MaxX = ceil(MaxZAxis/20) * 20;



Figure1=figure;
plot(zAxis,P,'LineWidth',2);
hold all;
plot(zAxis,Cl,'LineWidth',2);
plot(zAxis,Mg,'LineWidth',2);
plot(zAxis,K,'LineWidth',2);
hold off;
%title('Ion Distribution Profile','FontSize',TitleSize);
legend('P','Cl','Mg','K','Location','NorthWest');
xlabel('Z axis (angstrom)','FontSize',LabelSize);
ylabel('Concentration (M)','FontSize',LabelSize);
set(gca,'FontSize',AxisSize);

set(gca,'XLim',[MinX MaxX]);
set(gca,'XTick',[MinX:20:MaxX]);
%set(gca,'YLim',[0 9]);
%set(gca,'YTick',[0:3:9]);

saveas(Figure1,ionProfile,'epsc2')
