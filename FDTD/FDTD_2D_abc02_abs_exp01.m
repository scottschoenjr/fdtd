%**************************************************************************
%
% Finite difference time domain simulation (2D)
%
%   Two dimensional FDTD solution of the governing fluid equations. 
%
%
%            
% Change Log
%   Costas Arvanitis  2012      Initial version
%   Scott Schoen Jr   201612    Updates to excitation, plotting,
%                               organization
%
%**************************************************************************

close all;
clear all;
clc

% Visual settings
plotProgress = 1;
saveMovie = 1;
plotEvery = 2E-6; % Save a frame every __ seconds [s]

% 0 not to save, 1 to use default, or filename string
saveFilename = 1; 

% Load data with head properties
dataFile='..\data\head_2d_slc40_thr700_crop_res125.mat';
% dataFile='..\data\watercase2Dsmall_stratified.mat';
domain = load(dataFile);

%% Simulation parameters
f0 = 1E6;     % Pulse center frequency [Hz]
BW = 0.25;        % Pulse bandwidth
offset = 5E-6; % Pulse delay [s]

cw = 1479;         % Water sound speed [m/s]
rhow = 1000;       % Denstiy of Water [kg/m^3]
omega = 2.*pi.*f0; % Angular frequency [rad/s]

dx = domain.res*1E-3; % Distance between nodes [m]

% Receiver data
recPosition = 75E-3; % Distance from x = 0 to record data [m]
plotReceiverData = 1;

% Define sources
sourcePositions = [ ... % (x, y) positions [m]
    40, 30; 
    40, 31.5 ...
    ]./1E3;
sourceDelays = [0, 0]./1E6; % [s]

%% Load density, sound speed, and attneuation fields from data file
rho = flip( flip(domain.extrhotot, 2)', 2);   % Density
c = flip( flip(domain.extctot, 2)', 2);       % Sound speed
alpha = flip( flip(domain.extatttot, 2)', 2); % Attenuation

% Convert attenuation to diffusivity [m^2/sec]
delta = 0.0.*(f0/1E6)*(8/3)*alpha.*(c)./(omega^2);

xDim = size(rho, 2); % Number of nodes in x direction
yDim = size(rho, 1); % Number of nodes in y direction

xPositionVector = (0:dx:(xDim-1)*dx); % [m]
yPositionVector = (0:dx:(yDim-1)*dx); % [m]

% Initialize pressure and velocity fields
p = zeros(yDim, xDim, 3);
ux = zeros(yDim, xDim-1);
uy = zeros(yDim-1, xDim);

% Set Courrant compliant time step
dt = 0.55*dx/(sqrt(ndims(p))*max(c(:))*(1+omega^2*max(delta(:))^2/4));

% Pressure scalar difference coefficient
coP = (dt/dx)*rho.*c.^(2);

% Time to run the simulation
stopTime=1.8*xDim*dx/max(c(:));

% Create time vector
t = 0:dt:stopTime;

% If movie is to be saved, determine how many frames to skip between saved
% frames (e.g., if we have 25000 time steps, we probably don't need to save
% every single one).
if plotProgress
    
    % Create flag to determine if figure exsists
    progressPlotCreated = 0;
    
    % Get vectors for plotting receiver array position
    receiverArrayPlotX = [recPosition, recPosition ];
    receiverArrayPlotY = [ min(yPositionVector), max(yPositionVector) ];
    
    % Set movie parameters
    if saveMovie
        % Plot every __ frames
        plotEvery = max( 1, floor( plotEvery./dt ) ); 
        totalFrames= round(length(t)./plotEvery);
        % Initialize
        numFrames = 0;
        M(totalFrames) = struct('cdata',[],'colormap',[]);
    end
    
end

% Initialize receiver field
recPositionIndex = find( xPositionVector > recPosition, 1 );
aedata = zeros( yDim, length(t) );

% Position sources
numSources = length( sourcePositions(:, 1) );
excit_loc = 0.*sourcePositions;
for sourceCount = 1 : numSources
    
    % Define array of indices (rather than positions) of sources
    xCurrentSource = sourcePositions( sourceCount, 1 );
    yCurrentSource = sourcePositions( sourceCount, 2 );
    xIndex = find( xPositionVector > xCurrentSource,  1);
    yIndex = find( yPositionVector > yCurrentSource,  1);
    excit_loc( sourceCount, : ) = [xIndex, yIndex ];
    
end

% Define excitation pulse
pulseSignal = pulse(t, f0, BW, offset, 1);

%% Coefficients for 2nd order ABC
pOldup1 = zeros(3,xDim);
pOldup2 = zeros(3,xDim);
pOlddown1 = zeros(3,xDim);
pOlddown2 = zeros(3,xDim);
pOldright1 = zeros(yDim,3);
pOldright2 = zeros(yDim,3);
pOldleft1 = zeros(yDim,3);
pOldleft2 = zeros(yDim,3);

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

t(end)*1e+6
num=0.0;
ampl=1.0;
for n = 1:length(t)-2;
    
    % Assume the pulse ends early on, so we're not redefining the source
    % condition at every time step.
    if n <= length(t)/4
        for sourceCount = 1:numSources
            
            % Get delayed index
            delayedIndex = round( sourceDelays( sourceCount )./dt );
            
            % Get source location
            xIndexSource = excit_loc(sourceCount, 1);
            yIndexSource = excit_loc(sourceCount, 2);
            
            % Set pressure at the source (hard source)
            p(yIndexSource, xIndexSource, 1) = ...
                pulseSignal( n + delayedIndex );
            p(yIndexSource, xIndexSource, 2) = ...
                pulseSignal( n + 1 + delayedIndex );
            
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
          p(1:end, 2, 3) - pOldright2(:,2)...
        ) ...
        + abcCoefright(:,3).*pOldright1(:,2) - pOldright2(:,3);
    
    % ABC for left part of the domain
    p(1:end,end,3) = ...
        abcCoefright(:,1).*( p(1:end, end-2, 3) + pOldleft2(:,1) ) ...
        + abcCoefright(:,2).*( ...
          pOldleft1(:,1) + pOldleft1(:,3) - ...
          p(1:end, end-1, 3) - pOldleft2(:,2) ...
        )...
        + abcCoefright(:,3).*pOldleft1(:,2) - pOldleft2(:,3);
    
    % Update ABC parameters for next time step
    for kk=1:3
        
        pOldup2(kk, :) = pOldup1(kk,:);
        pOldup1(kk, :) = p(kk,1:end,3);
        
        pOlddown2(kk, :) = pOlddown1(kk,:);
        pOlddown1(kk, :) = p(end-kk+1, 1:end, 3);
        
        pOldright2(:, kk) = pOldright1(:, kk);
        pOldright1(:, kk) = p(1:end, kk, 3);
        
        pOldleft2(:, kk) = pOldleft1(:, kk);
        pOldleft1(:, kk) = p(1:end, end-kk+1, 3);
    end
    
    % Update stored pressure field
    p(:,:,1) = p(:,:,2);
    p(:,:,2) = p(:,:,3);
        
    % Plot field if desired
    if plotProgress && ( mod( n, plotEvery ) == 0 )
        
        % Create figure the first time
        if ~progressPlotCreated
            progressPlot = figure();
            hold all;
            progressPlotCreated = 1; % Set flag
        end

        % Plot data and format
        baselineData = ( c/max(max(c)) )./580;
        visPlotData = p(:,:,2) + baselineData;
        imagesc( 1E3.*xPositionVector, 1E3.*yPositionVector, ...
            visPlotData);
        
        % Plot receiver array position
        plot( receiverArrayPlotX.*1E3, receiverArrayPlotY.*1E3, ...
            '--w', 'LineWidth', 1.4 );
        
        axis image;
        xlabel('x Distance [mm]');
        ylabel('y Distance [mm]');
        titleString = sprintf( 'Time = %.2f us', n.*dt.*1E6);
        title(titleString, 'FontSize', 18, 'interpreter', 'tex');
        
        caxis([5E-5, 6E-2]);
        drawnow;
        
    end
    
    % Save data if we want to make a movie later
    if saveMovie
        if mod( n, plotEvery ) == 0
           numFrames = numFrames + 1;
           M( numFrames ) = getframe(gcf);
        end
    end
    
    % Save data at the receiver
    aedata(:, n)= p(:, recPositionIndex, 2)';
    
end;

% Plot received data if required
if plotReceiverData
    
    figure();
    
    % Create vectors and plot
    [tRecPlot, yRecPlot] = meshgrid( t, yPositionVector );
    pcolor( tRecPlot.*1E6, yRecPlot.*1E3, aedata );
    shading flat
    
    % Labels and formatting
    xlabel( 'Time [ms]' );
    ylabel( 'y Position [mm]' );
    
end

% Save the simulation data
if ~isequal( saveFilename, 0 )
   
    % Use specified filename, otherwise use default
    if isequal( saveFilename, 1 )
        currrentTime = clock;
        YYYY = num2str( clock(1) );
        MM = num2str( clock(2) );
        DD = num2str( clock(3) );
        hh = num2str( clock(4) );
        mm = num2str( clock(5) );
        ss = num2str( floor(clock(6)) );
        filenameString = [ ...
            YYYY, MM, DD, 'T', hh, mm, ss, '.mat' ...
            ];
    end
    
    % Save
    save( filenameString );
        
end
