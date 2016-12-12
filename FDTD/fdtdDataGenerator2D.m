%**************************************************************************
%
% Data Generator (FDTD)
%
%   Script to create data structure for use with FDTD code (2D). Data
%   structure is saved with variables
%   
%    res - Spacing of data points [mm]
%    zz - Vector of x positions [1-by-numX]
%    yy - Vector of y positions [1-by-numY]
%    extrhotot - Density at each (y, x) point [numY-by-numX]
%    extctot -  Sound speed at each (y, x) point [numY-by-numX]
%    extatttot - Attenuation at each (y, x) point  [numY-by-numX]
%
%              Scott Schoen Jr | Georgia Tech | 20161212
%
%************************************************************************** 

% Clear any existing variables
clear all;
clc;

% Set resolution !!in millimeters!!
res = 0.1250; % [mm]

% Set x- and z- vectors
xRange = [0, 200]./1E3; % [m]
yRange = [0, 80]./1E3; % [m]

% Get number of points in each direction
xx = xRange(1) : res./1E3 : xRange(2);
yy = yRange(1) : res./1E3 : yRange(2);

% Define density and sound speed profiles
[y, x] = meshgrid(yy, xx);
extrhotot = 1000.*ones( length(yy), length(xx) ); % [kg/m^3]
extctot = 1500.*ones( length(yy), length(xx) );   % [m/s]
extatttot = 1500.*ones( length(yy), length(xx) ); % [Np/m]

% Prompt user for filename
[saveFilename, saveDirectory ] = uiputfile( '*.mat', 'Save 2D data as:' );

% Save file
save( [saveDirectory, saveFilename] ); 

