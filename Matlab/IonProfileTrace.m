function IonProfileTrace(name)


TitleSize=16;
LabelSize=14;
AxisSize=12;
window=1000;



%load input data
ionFile1=sprintf('%s_P_profile_trace.dat', name);
File1=load(ionFile1);

ionFile2=sprintf('%s_Cl_profile_trace.dat', name);
File2=load(ionFile2);

ionFile3=sprintf('%s_Mg_profile_trace.dat', name);
File3=load(ionFile3);

ionFile4=sprintf('%s_K_profile_trace.dat', name);
File4=load(ionFile4);


ionProfileTrace=sprintf('%s_IonProfileTrace', name);



%extract variables
File1_in = File1';
File2_in = File2';
File3_in = File3';
File4_in = File4';

File1_in_cut = File1_in([1:5 end-4:end],:);
File2_in_cut = File2_in([1:5 end-4:end],:);
File3_in_cut = File3_in([1:5 end-4:end],:);
File4_in_cut = File4_in([1:5 end-4:end],:);

File1_bulk_trace = [];
for i=2:length(File1_in_cut)
    File1_bulk_trace = [File1_bulk_trace,mean(File1_in_cut(:,i))];
end
File1_bulk_trace_reduce = blockReduce(File1_bulk_trace,window);


File2_bulk_trace = [];
for i=2:length(File2_in_cut)
    File2_bulk_trace = [File2_bulk_trace,mean(File2_in_cut(:,i))];
end
File2_bulk_trace_reduce = blockReduce(File2_bulk_trace,window);


File3_bulk_trace = [];
for i=2:length(File3_in_cut)
    File3_bulk_trace = [File3_bulk_trace,mean(File3_in_cut(:,i))];
end
File3_bulk_trace_reduce = blockReduce(File3_bulk_trace,window);


File4_bulk_trace = [];
for i=2:length(File4_in_cut)
    File4_bulk_trace = [File4_bulk_trace,mean(File4_in_cut(:,i))];
end
File4_bulk_trace_reduce = blockReduce(File4_bulk_trace,window);


time = [];
for i=0:(length(File1_bulk_trace)-1)
    time = [time,0.0024*i];
end

time_reduce = [];
for i=0:(length(File1_bulk_trace_reduce)-1)
    time_reduce = [time_reduce,0.0024*i*window];
end


%find the common x and y range for plotting
MinTime = min(time);
MaxTime = max(time);
MinX = floor(MinTime/10) * 10;
MaxX = ceil(MaxTime/10) * 10;





% Plotting
Figure1=figure;
subplot(2,2,1);
hold all;
%p1=patchline(time,File1_bulk_trace,'lineWidth',1,'edgecolor',[204/256 229/256 255/256],'edgealpha',0.2);
p1=plot(time,File1_bulk_trace,'LineWidth',1,'Color',[255/256 240/256 245/256]);
plot(time_reduce,File1_bulk_trace_reduce,'LineWidth',2,'Color',[65/256 105/256 225/256]);
title('P','FontSize',TitleSize);
xlabel('Time (ns)','FontSize',LabelSize);
ylabel('Conc. (M)','FontSize',LabelSize);
uistack(p1,'bottom');
uistack(gca,'top');
set(gca,'XLim',[MinX MaxX],'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'FontSize',AxisSize,'layer','top');
hold off;


subplot(2,2,2);
hold all;
%p2=patchline(time,File2_bulk_trace,'lineWidth',1,'edgecolor',[204/256 229/256 255/256],'edgealpha',0.2);
p2=plot(time,File2_bulk_trace,'LineWidth',1,'Color',[255/256 240/256 245/256]);
plot(time_reduce,File2_bulk_trace_reduce,'LineWidth',2,'Color',[65/256 105/256 225/256]);
title('Cl','FontSize',TitleSize);
xlabel('Time (ns)','FontSize',LabelSize);
ylabel('Conc. (M)','FontSize',LabelSize);
uistack(p2,'bottom');
uistack(gca,'top');
set(gca,'XLim',[MinX MaxX],'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'FontSize',AxisSize,'layer','top');
hold off;


subplot(2,2,3);
hold all;
%p3=patchline(time,File3_bulk_trace,'lineWidth',1,'edgecolor',[204/256 229/256 255/256],'edgealpha',0.2);
p3=plot(time,File3_bulk_trace,'LineWidth',1,'Color',[255/256 240/256 245/256]);
plot(time_reduce,File3_bulk_trace_reduce,'LineWidth',2,'Color',[65/256 105/256 225/256]);
title('Mg','FontSize',TitleSize);
xlabel('Time (ns)','FontSize',LabelSize);
ylabel('Conc. (M)','FontSize',LabelSize);
uistack(p3,'bottom');
uistack(gca,'top');
set(gca,'XLim',[MinX MaxX],'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'FontSize',AxisSize,'layer','top');
hold off;


subplot(2,2,4);
hold all;
%p4=patchline(time,File4_bulk_trace,'lineWidth',1,'edgecolor',[204/256 229/256 255/256],'edgealpha',0.2);
p4=plot(time,File4_bulk_trace,'LineWidth',1,'Color',[255/256 240/256 245/256]);
plot(time_reduce,File4_bulk_trace_reduce,'LineWidth',2,'Color',[65/256 105/256 225/256]);
title('K','FontSize',TitleSize);
xlabel('Time (ns)','FontSize',LabelSize);
ylabel('Conc. (M)','FontSize',LabelSize);
uistack(p4,'bottom');
uistack(gca,'top');
set(gca,'XLim',[MinX MaxX],'XTick',[linspace(MinX,MaxX,6)]);
set(gca,'FontSize',AxisSize,'layer','top');
hold off;



saveas(Figure1,ionProfileTrace,'epsc2');
