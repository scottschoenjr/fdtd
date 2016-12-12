%**************************************************************************
%
% Finite difference time domain simulation (3D)
%
%   Code simulates the propagation of a microbubble emission (pulse) in the
%   brain. Tissue and skull dimensions are taken from CT scan data.
%
%
%
%            © Costas Arvanitis BWH, Harvard University 2012
%
%**************************************************************************

close all;
clear all;
clc

% Visual settings
plotProgress = 1;
saveMovie = 1;
plotEvery = 1E-6; % Save a frame every __ seconds [s]

% Load data with head properties
dataFile='..\data\head_2d_slc40_thr700_crop_res125.mat';
dataFile='..\data\watercase2D.mat';
domain = load(dataFile);

%% Some user modifiable parameters
son=14;
bps='n';                            %Choose exitation pulse "y' is from bubble dynamics 'n' is guasian pulse
% trgt=1;                             %Select the location of the point sourec
f0=1.32e6;       % Pulse center frequency [Hz]
BW = 0.4;        % Pulse bandwidth
offset = 2.5E-6; % Pulse delay [s]
cw = 1479;   % Water sound speed [m/s]
rhow = 1000; % Denstiy of Water [kg/m^3]
w = f0*2*pi;                        %Angular frequency
% lc = c0/f0;                       %Wavelength [m]
intconst = 8;                         %data points to collect images for video
dx = domain.res*1E-3; % Distance between nodes [m]


%% Initialise the pressure and velocity fields and set coefficients + a Courrant compliant dt
rho = flip( flip(domain.extrhotot, 2)', 2); % Density
c = flip( flip(domain.extctot, 2)', 2); % Sound speed
alpha = flip( flip(domain.extatttot, 2)', 2); % Attenuation
alpha(find(alpha==2.9))=2.9;%7.5x is correction to match the shear viscosity of brain tissue (3.5 kg-1sec-1)found in literature (the 2x is for stability)
alpha(find(alpha==0.000288))=0.000288;%30x is correction to mathc the shear viscosity of water (0.9e-4 kg-1sec-1)found in literature
% mi=alpha.*rho.*(c.^3)/w^2;% Shear viscosity
delta=0.0*(f0/1e+6)*(8/3)*alpha.*(c)./(w^2);% convert attenuation to diffusivity [m^2/sec]: 20 times smaller mi was used

xDim = size(rho, 2); % Number of nodes in x direction
yDim = size(rho, 1); % Number of nodes in y direction

% Initialize pressure and velocity fields
p = zeros(yDim, xDim, 3);
ux = zeros(yDim, xDim-1);
uy = zeros(yDim-1, xDim);

% Set Courrant compliant time step
dt = 0.55*dx/(sqrt(ndims(p))*max(c(:))*(1+w^2*max(delta(:))^2/4));

% Pressure scalar difference coefficient
coP = (dt/dx)*rho.*c.^(2);

% Time to run the simulation
stopTime=1.8*xDim*dx/max(c(:));

% Create time vector
t = 0:dt:stopTime;

% If movie is to be saved, determine how many frames to skip between saved
% frames (e.g., if we have 25000 time steps, we probably don't need to save
% every single one).
if saveMovie
    plotEvery = max( 1, floor( plotEvery./dt ) ); % Plot every __ frames
    totalFrames= round(length(t)./plotEvery);
    % Initialize
    numFrames = 0;
    M(totalFrames) = struct('cdata',[],'colormap',[]);
end

% Initialize receiver field
aedata = zeros( yDim, length(t) );

% Define excitation pulse
excit_loc = round( [xDim/2, yDim/2 ] );
numSources = length( excit_loc(:, 1) );
pulseSignal = pulse(t, f0, BW, offset, bps, 1);

%% Coefficients for 2nd order ABC
pOldup1=zeros(3,xDim);
pOldup2=zeros(3,xDim);
pOlddown1=zeros(3,xDim);
pOlddown2=zeros(3,xDim);
pOldright1=zeros(yDim,3);
pOldright2=zeros(yDim,3);
pOldleft1=zeros(yDim,3);
pOldleft2=zeros(yDim,3);

cr=2;

temp1=(dt/dx)*(c(1,:)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefup=zeros(3,xDim);
abcCoefup(1,:)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefup(2,:)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefup(3,:)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=(dt/dx)*(c(end,:)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefdown=zeros(3,xDim);
abcCoefdown(1,:)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefdown(2,:)=-2.0*(temp1-1.0./temp1)./temp2;
abcCoefdown(3,:)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=(dt/dx)*(c(:,end)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefright=zeros(yDim,3);
abcCoefright(:,1)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefright(:,2)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefright(:,3)= 4.0*(temp1+1.0./temp1)./temp2;

temp1=(dt/dx)*(c(:,1)*cr);
temp2=1.0./temp1+2.0+temp1;
abcCoefleft=zeros(yDim,3);
abcCoefleft(:,1)=-(1.0./temp1-2.0+temp1)./temp2;
abcCoefleft(:,2)=-2.0*(temp1 - 1.0./temp1)./temp2;
abcCoefleft(:,3)= 4.0*(temp1+1.0./temp1)./temp2;

timepnts=1:5*intconst:length(t)-1;
timestmps=zeros(length(timepnts),xDim);

t(end)*1e+6
num=0.0;
ampl=1.0;
for n = 1:length(t)-2;
    
    % Assume the pulse ends early on, so we're not redefining the source
    % condition at every time step.
    if n <= length(t)/4
        for sourceCount = 1:numSources
            
            % Get source location
            xIndexSource = excit_loc(sourceCount, 1);
            yIndexSource = excit_loc(sourceCount, 2);
            
            % Set pressure at the source (hard source)
            p(yIndexSource, xIndexSource, 1) = pulseSignal(n);
            p(yIndexSource, xIndexSource, 2) = pulseSignal(n+1);
            
        end
    end
    
    
    % Compute y-velocity
    uy = uy - ( 2*(dt/dx)./( rho(1:end-1,:)+rho(2:end,:) ) ).*(...
        ( 1 + ( (delta(1:end-1,:)+delta(2:end,:) )./(2*dt))).* ...
        ( p(2:end, :,2)- p(1:end-1,:,2) ) - ...
        ( (delta(1:end-1,:) + delta(2:end,:) )./(2*dt) ).*...
        ( p(2:end, :,1) - p(1:end-1,:,1)) ...
        );
    
    % Compute z-velocity
    ux = ux - ( 2*(dt/dx)./(rho(:,1:end-1)+rho(:,2:end) ) ).*...
        ( ( 1 + ( (delta(:,1:end-1) + delta(:,2:end) )./(2*dt))).* ...
        ( p(:,2:end,2) - p(:,1:end-1,2) ) - ...
        ( ( delta(:,1:end-1) + delta(:,2:end) )./(2*dt)).* ...
        ( p(:,2:end,1) - p(:,1:end-1,1) ) ...
        );
    
    % Update pressure field
    p(:,:,3) = p(:,:,2) - coP.*( ...
        [ uy; zeros(1, xDim) ] - [zeros(1,xDim); uy] + ...
        [ ux, zeros(yDim, 1) ] - [zeros(yDim,1), ux] ...
        );
    
    %% Define the 2nd order ABC
    
    % ABC for uper part of the domain
    p(1,1:end,3) = ...
        abcCoefup(1,:).*( p(3,1:end,3) + pOldup2(1,:) ) ...
        + abcCoefup(2,:).*( ...
        pOldup1(1,:) + pOldup1(3,:) - ...
        p(2, 1:end, 3) - pOldup2(2,:) ...
        ) ...
        + abcCoefup(3,:).*pOldup1(2,:) - pOldup2(3,:);
    
    % ABC for lower part of the domain
    p(end,1:end,3) = ...
        abcCoefdown(1,:).*( p(end-2,1:end,3) + pOlddown2(1,:) ) ...
        + abcCoefdown(2,:).*( ...
        pOlddown1(1,:) + pOlddown1(3,:) - ...
        p(end-1, 1:end, 3) - pOlddown2(2,:) ...
        )...
        + abcCoefdown(3,:).*pOlddown1(2,:) - pOlddown2(3,:);
    
    % ABC for right part of the domain
    p(1:end,1,3) = ...
        abcCoefright(:,1).*( p(1:end, 3, 3) + pOldright2(:,1) ) ...
        + abcCoefright(:,2).*( ...
        pOldright1(:,1) + pOldright1(:,3) - ...
        p(1:end, 2, 3) - pOldright2(:,2))   ...
        + abcCoefright(:,3).*pOldright1(:,2) - pOldright2(:,3);
    
    % ABC for left part of the domain
    p(1:end,end,3) = ...
        abcCoefright(:,1).*( p(1:end, end-2, 3) + pOldleft2(:,1) ) ...
        + abcCoefright(:,2).*( ...
        pOldleft1(:,1) + pOldleft1(:,3) - ...
        p(1:end, end-1, 3) - pOldleft2(:,2) ...
        )...
        + abcCoefright(:,3).*pOldleft1(:,2) - pOldleft2(:,3);
    
    % Update stored fields
    for kk=1:3
        
        pOldup2(kk,:) = pOldup1(kk,:);
        pOldup1(kk,:) = p(kk,1:end,3);
        
        pOlddown2(kk,:) = pOlddown1(kk,:);
        pOlddown1(kk,:) = p(end-kk+1,1:end,3);
        
        pOldright2(:,kk) = pOldright1(:,kk);
        pOldright1(:,kk) = p(1:end,kk,3);
        
        pOldleft2(:,kk) = pOldleft1(:,kk);
        pOldleft1(:,kk) = p(1:end,end-kk+1,3);
    end
    
    % Update stored pressure field
    p(:,:,1) = p(:,:,2);
    p(:,:,2) = p(:,:,3);
    
    
    % Plot field
    if plotProgress
        if n == 1
            figure();
        end
        
        xVisPlotVector = 0:dx:(xDim-1)*dx*1000;
        yVisPlotVector = 0:dx:(yDim-1)*dx*1000;
        visPlotData = flip( p(:,:,2) + c/max(max(c))/580, 2);
        imagesc(xVisPlotVector, yVisPlotVector, visPlotData);
        
        caxis([5E-5, 6E-2]);
        
        axis image;
        xlabel('Distance [mm]');
        set(gca, 'XDir', 'Reverse')
        ylabel('Distance [mm]')
        titleString = sprintf( 'Time = %.2f usec', n.*dt.*1E6);
        title(titleString, 'FontSize', 22);
        drawnow;
    end
    
    % Save data if we want to make a movie later
    if saveMovie
        if mod( n, plotEvery ) == 0
           numFrames = numFrames + 1;
           M( numFrames ) = getframe(gcf);
        end
    end
    
    %                 if isempty(find(n==1:5*intconst:length(t)-1, 1))
    %                     a=1;
    %                 else
    %                     currentframe = frame2im(getframe(fig)); % convert fig into image data
    %                     %     currentframe = im2bw(currentframe,0.4); % make it binary to save space
    %                     if n == 1
    %                         imwrite(currentframe, 'attscabs_1bps_trgt01_int02.tif', 'WriteMode', 'overwrite');
    %                     else
    %                         imwrite(currentframe, 'attscabs_1bps_trgt01_int02.tif', 'WriteMode', 'append');
    %                     end
    %                     close all;
    %                 end
    %
    % Save data at the reiver
    aedata(:,n)= p(:,end-5,2)';
    
end;

%% waterfall plots at the point source line coordinates
% figure
% axes('FontSize',16)
% simpleWaterfall(timestmps, 1, 1e4) % vertical offset = 1, scale factor = 1.9
% xlabel('Space [spatial index]','FontSize',16,'FontWeight','bold')
% ylabel('Time [frame number]','FontSize',16,'FontWeight','bold')
% axis([0 xDim -1 length(timepnts)]);
% figure
% Waterfall((1:xDim)*dx*1e3,t(timepnts)*1e6,timestmps) % vertical offset = 1, scale factor = 1.9
% colormap(gray)
% view(0,65)
% % axes('FontSize',16)
% xlabel('Space (mm)','FontSize',16,'FontWeight','bold')
% ylabel('Time (usec)','FontSize',16,'FontWeight','bold')
% zlabel('Amplitude','FontSize',16,'FontWeight','bold')
% grid off
% axis([0 (xDim*dx-dx)*1e3 t(1) t(end-1)*1e6 0 max(max(timestmps))]);

%% save the simulation data
% folder1=cd;
% folder1=Sdir;
% folder2=['1320_aeattscabs_2DHead_res' num2str(domain.res*1000) 'um_1gps0075_sontrgt0' num2str(son)]
% savefile =[folder1 '\' folder2 '.mat'];
% save(savefile,'aedata','dx','t','xDim','yDim','domain','excit_loc','timestmps','ratio','num','f0','son');