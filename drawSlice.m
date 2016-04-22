
function drawSlice(filename,Nx,Ny,Nz,angle,angle2)
% Takes in a PME or other grid and renders a couple of plots. Also creates 
% an average of the channel of a nanopore. Nx, Ny, and Nz each
% refer to the number of points along a dimension of the input grid; angle
% and angle2

% % READ CHANNEL GEOMETRY:
%channel=load('/home/alek/HEMOLYSIN/ANALYSIS/IONS/channelFluctuations.dat');

% 
% scaleZ = 147.959/Nz;
figure
scaleZ = 0.948769; % Taken directly from the dx file header

% zUp = channel(32,1);
% ZDown = channel(1,1);
% zStep = channel(2,1)-channel(1,1);
 
%zUp = 92;
%zDown = -14;
%zStep = 1;

v=zeros(Nx,Ny,Nz);

% Get the data from the dx file, disregarding the 8-line header.
[a, ~, ~] =importdata(filename, ' ', 8);

% Put the data in such a format that the next part will work.
b = a.data';  

tmp=1;
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            v(i,j,k)=b(tmp);
            tmp=tmp+1;
        end
    end
end


[x,y,z] = meshgrid(0:Ny-1,0:Nx-1,0:Nz-1);
    xmin = min(x(:)); 
    ymin = min(y(:)); 
    zmin = min(z(:));

    xmax = max(x(:)); 
    ymax = max(y(:)); 
    zmax = max(z(:));


% define planes here: 
    
centerX = Nx/2-2;
centerY = Ny/2-2;
centerZ = Nz/2-2;
d1 = cos(angle*pi/180);
d2 = sin(angle*pi/180);
d3 = 0.0;   
e1 = 0.0;
e2 = 0.0;
e3 = 1.0;
    
% n describes the resolution of the plot
n = 96;  
u = (0:0.5:n);
t = (-n:0.5:n);
umin = min(u(:)); 
tmin = min(t(:)); 

umax = max(u(:)); 
tmax = max(t(:)); 
[ii,Ni] = size(u);
[ii,Nj] = size(t);
for i=1:Ni
   for j=1:Nj
    X(i,j) = centerX + u(i)*d1 + t(j)*e1;
    Y(i,j) = centerY + u(i)*d2 + t(j)*e2;
    Z(i,j) = centerZ + u(i)*d3 + t(j)*e3;
   end 
end
hslice = surf(X,Y,Z);
    
xd1 = get(hslice,'XData');
yd1 = get(hslice,'YData');
zd1 = get(hslice,'ZData');
delete(hslice);



%figure('units','normalized','position',[.1 .1 .6 .6])
h = slice(x,y,z,v,xd1,yd1,zd1);

set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
hold on

%%%%%%%%%  plane II %%%%%%%%%

d1 = cos(angle2*pi/180);
d2 = sin(angle2*pi/180);
d3 = 0.0;   
e1 = 0.0;
e2 = 0.0;
e3 = 1.0;
    
    
n = 96;  
u = (0:0.5:n);
t = (-n:0.5:n);
umin = min(u(:)); 
tmin = min(t(:)); 

umax = max(u(:)); 
tmax = max(t(:)); 
[ii,Ni] = size(u);
[ii,Nj] = size(t);
for i=1:Ni
   for j=1:Nj
    X(i,j) = centerX + u(i)*d1 + t(j)*e1;
    Y(i,j) = centerY + u(i)*d2 + t(j)*e2;
    Z(i,j) = centerZ + u(i)*d3 + t(j)*e3;
   end 
end
hslice = surf(X,Y,Z);
    
xd2 = get(hslice,'XData');
yd2 = get(hslice,'YData');
zd2 = get(hslice,'ZData');
delete(hslice);


h2 = slice(x,y,z,v,xd2,yd2,zd2);
set(h2,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);

%h2x=contourslice(x,y,z,v,48,[],[],0.0);
%set(h3x,'EdgeColor',[1.0,1.0,1.0],'LineWidth',.5)
%%%%%%%%%

%hslice = surf(linspace(64,134,100),linspace(0,134,100),zeros(100));
 %     rotate(hslice,[-1,0,0],90.0, [0.0 48 48])
  %   rotate(hslice,[0,0,1],-45.0, [0.0 0.0 0.0]) 
%xd = get(hslice,'XData');
%yd = get(hslice,'YData');
%zd = get(hslice,'ZData');
%delete(hslice);
%
%h3 = slice(x,y,z,v,xd,yd,zd);
%set(h3,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8)



%hx = slice(x,y,z,v,48,[],[]);
%set(hx,'FaceColor','interp','EdgeColor','none')
h1x=contourslice(x,y,z,v,xd1,yd1,zd1,linspace(-0.5,0.3,25));
set(h1x,'EdgeColor',[0.0,0.0,0.0],'LineWidth',.5)

%daspect([1,1,1])
%axis off
%box on
%view(angle,0)
%camzoom(1.)
%camproj Orthographic

%%lightangle(-100,0)
%colormap (jet(20))
%set(gcf,'Renderer','zbuffer')

%colormap (flipud(jet(24)))
%caxis([-0.5,0.4])
%colorbar('vert')

h2x=contourslice(x,y,z,v,xd2,yd2,zd2,linspace(-0.5,0.3,25));
set(h2x,'EdgeColor',[0.0,0.0,0.0],'LineWidth',.5)

%contourslice(x,y,z,wind_speed,[xmin,100,xmax],ymax,zmin);
%set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)

%hy = slice(x,y,z,v,[],96,[]);
%set(hy,'FaceColor','interp','EdgeColor','none')

%hz = slice(x,y,z,v,[],[],96);
%set(hz,'FaceColor','interp','EdgeColor','none')

daspect([1,1,1])
axis off
box on
view(angle,0)
camzoom(1.)
camproj Orthographic

%lightangle(-100,0)
colormap (jet(20))
set(gcf,'Renderer','zbuffer')

colormap (flipud(jet(24)))
caxis([-0.5,0.3])
colorbar('vert')


%%%%%%% Getting average potential %%%%%%%%%


filePot = [filename '.potZ'];
fid = fopen(filePot,'w');
iCenter = Nx/2-1;
jCenter = Ny/2-1;

rOut = Nx/2 -2;
deltaTolerance =0.03; 
deltaR = 0.5;
i = iCenter;
j = jCenter;

% Cycle through the z-values and calculate the average value for the slice
% at that value. Treat one section differently as the inside of the
% nanopore.
for k=1:Nz 
      
    % This line specifically calibrated to start us at a given z.
    z = (k-Nz/2 + 1)*scaleZ + 29;  
           
    if  (( z > 82) || ( z < -26))    % True when above or below the channel
        average = 0.0;
        number = 0;
        for i=1:Nx
            for j=1:Ny
                if ( (i-iCenter)^2 + (j-jCenter)^2 < rOut^2 )
                    average = average + v(i,j,k);   
                    number = number +1; 
                end 
            end
        end
        
        fprintf(fid,'%g  %g %g \n', z, average/number, rOut);
%         
%     elseif ( z < -25)    % True when below the channel
%            
%         average = 0.0;  
%         number = 0;
%          for i=1:Nx
%             for j=1:Ny
%                 if ( (i-iCenter)^2 + (j-jCenter)^2 < (Nx/3)^2 )
%                     average = average + v(i,j,k);   
%                     number = number +1; 
%                 end 
%             end 
%         end 
%         fprintf(fid,'%g  %g %g\n', z, average/number, Nx/3.0);   
        
    else    % When inside the nanopore channel   
        average= v(iCenter,jCenter,k);     
        number =1; 
        
%         % Here I set the minimum R to a relatively high value to make the 
%         % cylindrical average the average of a hollow cylinder.

        % Here I set the minimum R to a low value to preven the while loop
        % from terminating before taking a good average. 
        R=0.75; deltaAverage = 0.0;
     
 %       while (R < 9)
        while (abs(deltaAverage) < deltaTolerance) & (R < Nx/2-1)
            % Compute a cylindrical average of the channel, 
            % increasing the radius of averaging each step. If the average
            % changes too dramatically on any step, quit and put out the
            % old average.
            R = R+deltaR;
            averageOld = average;
            numberOld = number;
            average = 0.0;
            number = 0;
            for i=1:Nx
                for j=1:Ny
                    if ( (i-iCenter)^2 + (j-jCenter)^2 < R^2 )
                        average = average + v(i,j,k);   
                        number = number +1; 
                    end 
                end 
            end 
            average = average/number;        
        
            if (number ~= numberOld)  % Avoid divide by 0 errors
                
                % Gives the fractional difference between the new and the old averages
                deltaAverage = (average*number -averageOld*numberOld)/(number-numberOld)...
                -averageOld;   
                
                % Gives a numeric difference between the new and the old averages
                %deltaAverage = average-averageOld  
            end 
        end
        fprintf(fid,'%g  %g %g \n', z, averageOld, R-deltaR);              
     end
      
end
status = fclose(fid);
