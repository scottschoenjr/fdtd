%**************************************************************************
%
% Create Localized Excitation
%
%   Function creates a Gaussian-pulse excitation with input parameters for
%   its initial time, central frequency, bandwidth, and time step.
%
% Inputs
%   t         - Time vector [s]
%   dt        - Time Step [s]
%   f0        - Center frequency of pulse [Hz]
%   BW        - Fractional bandwidth of pulse
%   offset    - Time offset from center of signal
%   plotPulse - Display waveform and spectrum (plotted if == 1)
%
% Returns
%   pulseSignal - Time series pulse
%
% Change Log
%   2013XXXX Dr Costas Arvanitis   Original
%   201610XX Scott Schoen Jr       Formatting, added plotting option
%
%**************************************************************************

function pulseSignal = pulse( timeVector, f0, BW, offset, plotPulse)

% Apply shift
t = timeVector - offset;

% If not specified, do not plot pulse
if nargin < 5
    plotPulse = 0;
end

% Define Gaussian pulse excitation here. This function replaces the 
% gauspuls function to minimize required toolboxes. We'll define parameters
% like the BW ratio, even though they're not adjustable from the function
% call.
bwRatio = -6;  % [dB] Determines how quickly envelope falls off
freqVariance = ...
    -( BW.^(2).*f0.^(2) )./ ( 8.*(bwRatio./20) ); % Freq. variance [s^-2]
timeVariance = 1/(4.*pi.^(2).*freqVariance);      % Time Variance [s^2]
pulseEnvelope = exp(-t.^(2)/(2.*timeVariance));
pulseSignal = pulseEnvelope.*cos( 2.*pi.*f0.*t );

% Ensure maximum is at +1
[maxValue, maxIndex] = max( abs(pulseSignal) );
if pulseSignal(maxIndex) < 0
    pulseSignal = -1.*pulseSignal;
end

% Plot pulse if prompted
if plotPulse
    
    % Determine what time scale to plot over
    numSigmas = 8;
    timeSigma = sqrt( timeVariance ); % [s]
    tStart = offset - numSigmas.*timeSigma;
    tEnd = offset + numSigmas.*timeSigma;
    
    % Get appropriate indices
    startIndex = find( timeVector > tStart, 1 );
    endIndex = find( timeVector > tEnd, 1 );
    if tStart < timeVector(1)
        startIndex = 1;
    end
    if tEnd > timeVector(end)
        endIndex = length(t);
    end
    
    
    figure()
    
    % Plot the time series
    subplot( 2, 1, 1 )
    plot( timeVector(startIndex:endIndex).*1E3,  ...
        pulseSignal(startIndex:endIndex), 'k' );
    xlabel( 'Time [ms]' );
    ylabel( 'Amplitude [AU]' );
    
    % Plot the frequency content
    subplot( 2, 1, 2 )
    dt = t(2) - t(1); % Sample step [s]
    Fs = 1./dt;       % Sampling frequency [Hz]
    fVector = linspace( 0, Fs, length(t(startIndex:endIndex)) );
    spectrum = abs( fft( pulseSignal(startIndex:endIndex) ) );
    spectrumNorm = spectrum./max(spectrum);
    plot( fVector./1E3, spectrumNorm, 'k' );
    xlabel( 'Frequency [kHz]' );
    xlim( [0, (Fs./2)./1E3] );
    ylabel( 'Normalized Amplitude [AU]' );
    
end

end

