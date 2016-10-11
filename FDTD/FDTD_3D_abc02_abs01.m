%**************************************************************************
%
% Finite difference time domain simulation
%
%   Code simulates the propagation of a microbubble emission (pulse) in the
%   brain. Tissue and skull dimensions are taken from CT scan data.
%
%
%
%             © Costas Arvanitis BWH, Harvard University 2013 
%
%**************************************************************************

clear all;
close all;
clc;

% Inputs ------------------------------------------------------------------
bps='n';        % Choose exitation pulse "y' is from bubble dynamics
trgt = 5;       % select the location of the point sourec
f0 = 0.88e6;    % Mean frequency of the excitation pulse in [Hz]
w = f0*2*pi;    % Angular frequency [rad/s]
BW = 0.6;       % Fractional bandwidth of source 
cw = 1479;      % Speed of sound of water at 20 Celcius [m/s]
rhow = 1000;    % Denstiy of Water [kg/m^3]

intconst = 1;    % Data points to collect images for video
% -------------------------------------------------------------------------

% Plot options
plotSkull = 1;
waterSimulation = 0;
plotPropagation = 1;

% Add necessary paths
addpath( genpath( '..\pam_recon' ) );

% Load density, sound speed, and attenutaion parameters from data file
dataFile = ...
    '..\..\orig\data\simulated\head_3dproperties_thr700_extcrop04_int01.mat';
domain = load( dataFile );

%% Initialise the p and u fields and set coefficients + a Courrant compliant dt
rho = domain.rhotot; % Density profile [kg/m^3]
dx = domain.dvi*1E-3; % Length of nodes in x direction [m]

% For Sim with skull and brain
c = domain.ctot;
alpha = domain.atttot;
alpha(alpha==2.9)=2.9;%7.5x (1/3)is correction to match the shear viscosity of brain tissue (3.5 kg-1sec-1)found in literature
alpha(alpha==0.000288)=0.000288;%30x (1/3)is correction to match the shear viscosity of water (0.9e-4 kg-1sec-1) found in literature
% mi=alpha.*rho.*(c.^3)/w^2;% Shear viscosity
delta = (f0/1e+6)*(8/3)*alpha.*(c)./(w^2); % Convert attenuation to diffusivity [m^2/s]

% Get the size of the grid
xDim = size(rho,1);   % Number of nodes in x direction
yDim = size(rho,2);   % Number of nodes in y direction
zDim = size(rho,3);   % Number of nodes in z direction

p = zeros(xDim,yDim,zDim,3);          %pressure scalar field initialisation
ux = zeros(xDim,yDim-1,zDim);       %velocity vector field initialisation
uy = zeros(xDim-1,yDim,zDim);       %velocity vector field initialisation
uz = zeros(xDim,yDim,zDim-1);       %velocity vector field initialisation

% Set max simulation time (can be reduced for pure water case)
stopTime = 6.3*xDim*dx/max(c(:));

% For a simulation in water
if waterSimulation
    % Time to run the simulation
    stopTime = 3.34*xDim*dx/max(c(:));  
    % Create uniform density and sound speed fields
    c = cw*ones(xDim,yDim,zDim);
    rho= rhow*ones(xDim,yDim,zDim);
    % Set uniform terms
    delta =(f0/1e+6)*3.716e-14*ones(xDim,yDim,zDim);%*3.72e-14 or 1.486e-13
    ratio = max(c(:))/(1.5*f0)/(dx);
end

%% Display Skull if Desired
if ( plotSkull && ~waterSimulation )
    figure(101)
    xslice = [xDim/2, xDim]; 
    yslice = [yDim/2, yDim]; 
    zslice = [1,zDim/2];
    skullPlotHandle = slice(dx*(1:yDim), dx*(1:xDim), dx*(1:zDim), rho, ...
        dx*yslice, dx*xslice, dx*zslice);
    set(skullPlotHandle,'edgecolor','none')
    axis image;
    view([-70 10]);
    colormap gray(100)
    xlabel('Distance [mm]','FontWeight','bold')
    ylabel('Distance [mm]','FontWeight','bold')
    zlabel('Distance [mm]','FontWeight','bold')
end

% Set time step small enough to satisfy courrant condition
% dt = 0.55*dx/(sqrt(ndims(p))*max(c(:))*(1+w^2*max(delta(:))^2/4));
% Set time step small enough to satisfy courrant condition
dt = 0.55*dx/(sqrt(ndims(p))*max(c(:))); 

% Pressure scalar difference co-efficient
coP = (dt/dx)*rho.*c.^2;     

% Set up time vector
t = 0:dt:stopTime;

%% Define the excitation pulse
% excit_loc=domain.targets/(1e+3*dx); %location of multiple sources locations
c0 = mean(mean(mean(c)));
lambda = c0./f0;
deltaR = lambda; % Difference in space [m]
deltaPx = deltaR./dx; % Differece in pixels
loc1 = round( [xDim/2, yDim/2, zDim/2 ] );
loc2 = round( [xDim/2 - 2.*deltaPx, yDim/2, zDim/2 ] );
loc3 = round( [xDim/2 + 2.*deltaPx, yDim/2, zDim/2 ] );
loc4 = round( [xDim/2, yDim/2 - 2.*deltaPx, zDim/2 ] );
loc5 = round( [xDim/2, yDim/2 + 2.*deltaPx, zDim/2 ] );
excit_loc = round( [ loc1; loc2; loc3; loc4; loc5 ] ); %multiple sources locations
offset = 2.3;
excit_b = (1480*1000)/(1300*1562)*0.2*pulse(t, f0, BW, offset, bps, 0);  % 21.45   with bubbles  b7.868x with bubble excitation and 5.655 at 440,  4 at 880 kHz and  3.268 at 13200 kHz pulse
sum(abs(fft(excit_b))/length(t))
% fig=figure;
pdata=zeros(xDim,zDim,length(t));

%% Initilize vectors for coefficients for 2nd order Absorption Boundary Conditions
pOldup1      = zeros(3,yDim,zDim);
pOldup2      = zeros(3,yDim,zDim);
pOlddown1    = zeros(3,yDim,zDim);
pOlddown2    = zeros(3,yDim,zDim);
pOldright1   = zeros(xDim,3,zDim);
pOldright2   = zeros(xDim,3,zDim);
pOldleft1    = zeros(xDim,3,zDim);
pOldleft2    = zeros(xDim,3,zDim);
pOldforward1 = zeros(xDim,yDim,3);
pOldforward2 = zeros(xDim,yDim,3);
pOldback1    = zeros(xDim,yDim,3);
pOldback2    = zeros(xDim,yDim,3);

cr=2;

temp1 = squeeze((dt/dx)*(c(1,:,:)*cr));
temp2 = 1.0./temp1+2.0+temp1;
abcCoefup = zeros(3,yDim,zDim);
abcCoefup(1,:,:) = -(1.0./temp1-2.0+temp1)./temp2;
abcCoefup(2,:,:) = -2.0*(temp1- 1.0./temp1)./temp2;
abcCoefup(3,:,:) = 4.0*(temp1+1.0./temp1)./temp2;

temp1 = squeeze((dt/dx)*(c(end,:,:)*cr));
temp2 = 1.0./temp1+2.0+temp1;
abcCoefdown = zeros(3,yDim,zDim);
abcCoefdown(1,:,:) = -(1.0./temp1-2.0+temp1)./temp2;
abcCoefdown(2,:,:) = -2.0*(temp1-1.0./temp1)./temp2;
abcCoefdown(3,:,:) = 4.0*(temp1+1.0./temp1)./temp2;

temp1=squeeze((dt/dx)*(c(:,end,:)*cr));
temp2=1.0./temp1+2.0+temp1;
abcCoefright=zeros(xDim,3,zDim);
abcCoefright(:,1,:)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefright(:,2,:)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefright(:,3,:)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=squeeze((dt/dx)*(c(:,1,:)*cr));
temp2=1.0./temp1+2.0+temp1;
abcCoefleft=zeros(xDim,3,zDim);
abcCoefleft(:,1,:)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefleft(:,2,:)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefleft(:,3,:)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=squeeze((dt/dx)*(c(:,:,end)*cr));
temp2=1.0./temp1+2.0+temp1;
abcCoefforward=zeros(xDim,yDim,3);
abcCoefforward(:,:,1)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefforward(:,:,2)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefforward(:,:,3)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=squeeze((dt/dx)*(c(:,:,1)*cr));
temp2=1.0./temp1+2.0+temp1;
abcCoefback=zeros(xDim,yDim,3);
abcCoefback(:,:,1)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefback(:,:,2)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefback(:,:,3)= 4.0*(temp1+1.0./temp1)./temp2;

timepnts=1:5*intconst:length(t)-1;
t(end)*1e+6
for n = 1:length(t)-1;
    
    
    if n<=length(t)/3
        for sourceCount = 1:trgt
            
            % Determine current source location
            xIndexSource = excit_loc(sourceCount, 1, 1);
            yIndexSource = excit_loc(sourceCount, 2, 1);
            zIndexSource = excit_loc(sourceCount, 3, 1);
            
            % Set this source's excitation values at this and the 
            % next time step
            p( xIndexSource, yIndexSource, zIndexSource, 1) = ...
                excit_b(n);
            p( xIndexSource, yIndexSource, zIndexSource, 2) = ...
                excit_b(n + 1);
        end
    end
    
    % Compute the y velocity
    uy = uy - ( 2*(dt/dx)./(rho(1:end-1,:,:) + rho(2:end,:,:)) ).*...
        ( ...
          ( 1 + ( (delta(1:end-1,:,:) + delta(2:end,:,:) )./(2*dt)) ).*...
            ( p(2:end,:,:,2) - p(1:end-1,:,:,2) ) - ...
          ( (delta(1:end-1,:,:) + delta(2:end,:,:) )./(2*dt) ).*...
            ( p(2:end,:,:,1) - p(1:end-1,:,:,1) ) ...
         ); 
    % Compute the x velocity
    ux = ux - ( 2*(dt/dx)./(rho(:,1:end-1,:) + rho(:,2:end,:)) ).* ...
        ( ...
          ( 1 + ( (delta(:,1:end-1,:) + delta(:,2:end,:) )./(2*dt)) ).*...
            ( p(:,2:end,:,2)- p(:,1:end-1,:,2) ) - ...
          ( (delta(:,1:end-1,:) + delta(:,2:end,:) )./(2*dt) ).* ...
            ( p(:,2:end,:,1) - p(:,1:end-1,:,1) ) ...
          );
      
    % Compute the z-velocity
    uz = uz - ( 2*(dt/dx)./(rho(:,:,1:end-1) + rho(:,:,2:end)) ).*...
        ( ...
          ( 1 + ( (delta(:,:,1:end-1) + delta(:,:,2:end))./(2*dt)) ).*...
            p(:,:,2:end,2) - p(:,:,1:end-1,2) ) - ...
          ( (delta(:,:,1:end-1) + delta(:,:,2:end))./(2*dt)).* ...
            ( p(:,:,2:end,1)- p(:,:,1:end-1,1) ) ...
          ); 
    
    % Calculate pressure from velocity
    p(:,:,:,3) = p(:,:,:,2) - coP.*(...
          cat(1, uy, zeros(1, yDim, zDim)) - cat(1,zeros(1, yDim, zDim), uy) ...
        + cat(2, ux, zeros(xDim, 1, zDim)) - cat(2,zeros(xDim, 1, zDim), ux) ...
        + cat(3, uz, zeros(xDim, yDim, 1)) - cat(3,zeros(xDim, yDim, 1), uz) ...
        );   
    
    
    %% Apply 2nd order ABC
    
    % For the upper part of the domain...
    p(1, 1:end, 1:end, 3) = ...
          abcCoefup(1,:,:).*( p(3,1:end,1:end,3) + pOldup2(1,:,:) ) ...
        + abcCoefup(2,:,:).*( ...
              pOldup1(1,:,:) + pOldup1(3,:,:) - ...
              p(2,1:end,1:end,3) - pOldup2(2,:,:) ...
              )...
        + abcCoefup(3,:,:).*pOldup1(2,:,:) ...
        - pOldup2(3,:,:);
    
    % ...lower part...
    p(end, 1:end, 1:end, 3) = ...
          abcCoefdown(1,:,:).*( p(end-2,1:end,1:end,3) + pOlddown2(1,:,:) )...
        + abcCoefdown(2,:,:).*( ...
              pOlddown1(1,:,:) + pOlddown1(3,:,:) - ...
              p(end-1,1:end,1:end,3) - pOlddown2(2,:,:) ...
              )...
        + abcCoefdown(3,:,:).*pOlddown1(2,:,:) ...
        - pOlddown2(3,:,:);
    
    % ...right part...
    p(1:end,1,1:end,3) = ...
          abcCoefright(:,1,:).*( p(1:end,3,1:end,3) + pOldright2(:,1,:) )...
        + abcCoefright(:,2,:).*( ...
              pOldright1(:,1,:) + pOldright1(:,3,:) - ...
              p(1:end,2,1:end,3) - pOldright2(:,2,:) ...
              ) ...
        + abcCoefright(:,3,:).*pOldright1(:,2,:) ...
        - pOldright2(:,3,:);
    
    % ...left part...
    p(1:end, end, 1:end, 3) = ...
          abcCoefright(:,1,:).*( p(1:end,end-2,1:end,3) + pOldleft2(:,1,:) )...
        + abcCoefright(:,2,:).*( ...
              pOldleft1(:,1,:) + pOldleft1(:,3,:) - ...
              p(1:end,end-1,1:end,3) - pOldleft2(:,2,:) ...
              ) ...
        + abcCoefright(:,3,:).*pOldleft1(:,2,:) ...
        - pOldleft2(:,3,:);
    
    % ...forward part...
    p(1:end, 1:end, 1, 3) = ...
          abcCoefforward(:,:,1).*( p(1:end,1:end,3,3) + pOldforward2(:,:,1) ) ...
        + abcCoefforward(:,:,2).*( ...
              pOldforward1(:,:,1) + pOldforward1(:,:,3) - ...
              p(1:end,1:end,2,3) - pOldforward2(:,:,2) ...
              )...
        + abcCoefforward(:,:,3).*pOldforward1(:,:,2) ...
        - pOldforward2(:,:,3);
    
    % ..and back part of the domain.
    p(1:end,1:end,end,3) = ...
          abcCoefback(:,:,1).*( p(1:end,1:end,end-2,3) + pOldback2(:,:,1) ) ...
        + abcCoefback(:,:,2).*( ...
              pOldback1(:,:,1) + pOldback1(:,:,3) - ...
              p(1:end,1:end,end-1,3) - pOldback2(:,:,2) ...
              )...
        + abcCoefback(:,:,3).*pOldback1(:,:,2) ...
        - pOldback2(:,:,3);
    
    % Update stored fields
    for kk=1:3
        
        pOldup2(kk,:,:)=pOldup1(kk,:,:);
        pOldup1(kk,:,:)=p(kk,1:end,1:end,3);
        
        pOlddown2(kk,:,:)=pOlddown1(kk,:,:);
        pOlddown1(kk,:,:)=p(end-kk+1,1:end,1:end,3);
        
        pOldright2(:,kk,:)=pOldright1(:,kk,:);
        pOldright1(:,kk,:)=p(1:end,kk,1:end,3);
        
        pOldleft2(:,kk,:)=pOldleft1(:,kk,:);
        pOldleft1(:,kk,:)=p(1:end,end-kk+1,1:end,3);
        
        pOldforward2(:,:,kk)=pOldforward1(:,:,kk);
        pOldforward1(:,:,kk)=p(1:end,1:end,kk,3);
        
        pOldback2(:,:,kk)=pOldback1(:,:,kk);
        pOldback1(:,:,kk)=p(1:end,1:end,end-kk+1,3);
    end
    
    p(:,:,:,1) = p(:,:,:,2);
    p(:,:,:,2) = p(:,:,:,3);
    
    % Capture the data in the imaging plane 
    pdata(:,:,n) = squeeze( p(:, end-5, :, 3) );
    %% 3D data visualizatioon
    if plotPropagation
        v = p(:,:,:,3)*8e6+c;%+0.0*domain.mrvol;%data scaling
        xslice = [double(excit_loc(trgt,1,1)),xDim]; yslice =[1,double(excit_loc(trgt,2,1)),yDim]; zslice =[1,double(excit_loc(trgt,3,1)),zDim];
        hsrf=slice(1e3*dx*(1:yDim),1e3*dx*(1:xDim),1e3*dx*(1:zDim),v,1e3*dx*yslice,1e3*dx*xslice,1e3*dx*zslice);
        set(hsrf,'edgecolor','none','FaceAlpha',0.75)
        caxis([1e3 3e3])
        view([25 30]);
        colormap jet(100)
        axis([1e3*dx 1e3*dx*xDim 1e3*dx 1e3*dx*yDim 1e3*dx 1e3*dx*zDim])
        axis equal
        xlabel('Distances (mm)','FontWeight','bold')
        %             ylabel('Distance (mm)','FontWeight','bold')
        zlabel('Distance (mm)','FontWeight','bold')
        title('3D FDTD Simulation','FontSize',10,'FontWeight','bold')
        text(55,110,sprintf('Time = %.2f usec',n*dt*1e6),'FontSize',12,'FontWeight','bold','HorizontalAlignment','right')
        text(85,30,-12,'Costas Arvanitis','FontSize',10,'FontWeight','bold','HorizontalAlignment','right')
        text(100,25,-12,'Harvard Medical School, BWH','FontSize',10,'FontWeight','bold','HorizontalAlignment','right'); drawnow;
    end
    
    %% 2D visualisation
    %     subplot (1,3,1)
    %     imagesc(1e3*dx*(1:yDim),1e3*dx*(1:xDim),squeeze(v(:,:,double(excit_loc(4,3,1)))))
    %     set(gca,'YDir','Reverse')
    %     colormap(gray)
    %     caxis([500 4e3])
    %     axis image
    %     xlabel('Distance (mm)')
    %     ylabel('Distance (mm)')
    %     title('Coronal Plane','Color',[0 0 0],'FontSize', 12)
    %     subplot (1,3,2)
    %     imagesc(1e3*dx*(1:zDim),1e3*dx*(1:xDim),squeeze(flipdim(v(:,double(excit_loc(4,2,1)),:),1)))
    %     set(gca,'xDir','normal')
    %     colormap(gray)
    %     axis image
    %     caxis([500 4e3])
    %     xlabel('Distance (mm)')
    %     title('Sagittal Plane','Color',[0 0 0],'FontSize', 12)
    %     subplot (1,3,3)
    %     imagesc(1e3*dx*(1:yDim),1e3*dx*(1:zDim),squeeze(v(double(excit_loc(4,1,1)),:,:))')
    %     set(gca,'YDir','normal')
    %     caxis([500 4e3])
    %     colormap(gray)
    %     axis image
    %     xlabel('Distance (mm)')
    %     title('Transverse Plane','Color',[0 0 0],'FontSize', 12)
    %     text(50,110,sprintf('Time = %.2f usec',n*dt*1e6),'FontSize',12,'HorizontalAlignment','right');
    %     text(50,-30,'Video courtesy of Costas Arvanitis','FontSize',7,'HorizontalAlignment','right');
    %     text(50,-35,'BWH, Harvard University','FontSize',7,'HorizontalAlignment','right'); drawnow;
    %% Convert imges to tiff for making video
    %             if isempty(find(n==1:5*intconst:length(t)-1, 1))
    %                 a=1;
    %             else
    %                 currentframe = frame2im(getframe(fig)); % convert fig into image data
    %                 %     currentframe = im2bw(currentframe,0.4); % make it binary to save space
    %                 if n == 1
    %                     imwrite(currentframe, 'attscabs_1bps_trgt01_int02.tif', 'WriteMode', 'overwrite');
    %                 else
    %                     imwrite(currentframe, 'attscabs_1bps_trgt01_int02.tif', 'WriteMode', 'append');
    %                 end
    %                 close all;
    %             end
    %
    fprintf('\b\b\b%3.0f',100*n/length(t));
    
end;
aedata=pdata(:,:,(round(length(t)/4):round(length(t))));

%% display data
% passive detector rf data
% figure (3)
% [x,y,z] = meshgrid(1:zDim,1:xDim,round(length(t)/4):length(t));
% v = aedata;%pdata(:,:,round(length(t)/3):end);
% xslice = [1,double(excit_loc(2,3))]; yslice =[double(excit_loc(2,1)),xDim]; zslice =[length(t)/1.2,length(t)/2];
% h=slice(1e3*dx*x,1e3*dx*y,1e6*dt*z,v,round(1e3*dx*xslice),1e3*dx*yslice,1e6*dt*zslice);
% set(h,'edgecolor','none','FaceAlpha',0.75)
% axis([1e3*dx 1e3*dx*zDim 1e3*dx 1e3*dx*xDim 1e6*dt*length(t)/2 1e6*dt*length(t)])
% view([45 55]);
% colormap jet(100)
% xlabel('Distance (mm)','FontWeight','bold')
% ylabel('Distance (mm)','FontWeight','bold')
% zlabel('Time (usec)','FontWeight','bold')
% title('Passive Emissions Detector','Color',[0 0 0],'FontSize', 22); drawnow;

% %% save the simulation data
% folder1=Sdir;
% % folder1=cd;
% folder2=['880kHz_aeattscabs1_3Dwater04_int0' num2str(domain.interp) '_1gps_trgt0' num2str(trgt) 'asa'];
% savefile =[folder1 '\' folder2 '.mat'];
% dumvar=0;
% save(savefile,'dumvar','-v7.3');
% save(savefile,'aedata','dx','t','xDim','yDim','zDim','excit_loc','ratio','trgt','-v7.3');