%**************************************************************************
%
% Create Localized Excitation
%
%   Function creates a Gaussian-pulse excitation with input parameters for
%   its initial time, central frequency, bandwidth, and time step.
%
%
%
%

function excit_b = pulse( t, dt, f0, BW, offset, bps, plotPulse)

% If not specified, do not plot pulse
if nargin < 7
    plotPulse = 0;
end

% Define Gaussian pulse excitation
excit_b = gauspuls(t-t(end)/2, f0, BW);

% Compute the frequency vector
Fs = 1./dt;
numPoints = length(t);
startValue = -(numPoints./2) + 1 - 0.5;
endValue = length(t)/2 - 0.5;
omegaVector = (Fs).*2.*pi.*(startValue : endValue)./numPoints;

% Define the excitation vector, which is shifted to fit the simulation
% timescale
timeShift = (t(end))/4*offset;
phaseShift = -1j.*omegaVector.*timeShift;
% Get the shifted signal in the frequency domain
pulseFD = fftshift( fft(excit_b) ).*exp( phaseShift );
% Transform back to time domain and take real part
pulseTD = real( ifft( ifftshift( pulseFD ) ) );

% Return the excitation vector
excit_b = pulseTD;

[yy tt]=max(abs(excit_b));
if excit_b(tt)<0
aa=-1.0;
else
    aa=1.0;
end
excit_b=aa*excit_b;


% Plot pulse if prompted
if plotPulse
    
    % Determine what time scale to plot over
    tWidth = (1./(BW.*f0));
    tMiddle = ( t(end) - t(1) )./2;
    tStart = tMiddle + timeShift./1E3 - tWidth;
    tEnd = tMiddle + timeShift./1E3 + tWidth;
%     startIndex = find( t > tStart, t );
%     endIndex = find( t > tEnd, 1 );
    
    startIndex = 1;
    endIndex = length(t);
    
    
    figure()
    % Plot the time series
    subplot( 3, 1, 1 )
    plot( t(startIndex:endIndex).*1E3,  excit_b(startIndex:endIndex) );
    xlabel( 'Time [ms]' );
    ylabel( 'Amplitude [AU]' );
    
    % figure (1)
    % subplot(2,2,1);title('Time Trace','FontWeight','bold')
    % plot(1e6*t(1:600),0.2*excit_b(1:600),'k')
    % xlabel('Time (Sec)','FontWeight','bold')
    % ylabel('Amplitude (a.u.)','FontWeight','bold')
    % hold on
    % subplot(2,2,2);title('Power Spectrum','FontWeight','bold')
    % freq=(0:length(t))/dt/length(t);
    % ampl=0.2*(abs(fft(excit_b))/length(t));
    % plot(1e-6*freq(1:500),ampl(1:500),'k')
    % xlabel('Frequency (kHz)','FontWeight','bold')
    % ylabel('Amplitude (a.u.)','FontWeight','bold')
    % sum(ampl)
    % hold on
    
end
if (bps=='y')         
bbl=load('C:\work_material\data\in_vivo\monkey\PUS\06232012\nhp01\Registered_data\mr_ct_simdata\analysis\point sources\BubbleSim_pulse_440.mat');
dt1=bbl.simulation.tr(2);
excit_b=bbl.simulation.pr(1:(dt/dt1):end);
excit_b=[excit_b;excit_b(end)*ones(length(t)-length(excit_b),1)];
end
% excit_b=excit_b/sum(abs(excit_b));

% figure
% subplot(2,2,3);title('Time Trace','FontWeight','bold')
% plot(1e6*t(1:600),4.2364*excit_b(1:600),'k')
% xlabel('Time (usec)','FontWeight','bold')
% ylabel('Presure (Pa)','FontWeight','bold')
% subplot(2,2,4);title('Power Spectrum','FontWeight','bold')
% freq=(0:length(t))/dt/length(t);
% ampl=4.2364*(abs(fft(excit_b))/length(t));
% plot(1e-6*freq(1:500),ampl(1:500),'k')
% xlabel('Frequency (kHz)','FontWeight','bold')
% ylabel('Amplitude','FontWeight','bold')
% sum(ampl)
% Atempt to create sinusoidal pulse
% A=2;
% f= 500000;
% phi=0;
% a=100000;
% 
% dt=.0000001;
% f=f*2*pi;
% t=0:dt:dt*1e3;
% y1=A*sin(f*t + phi).*exp(-a*t)+A/2*sin(2*f*t + phi).*exp(-a*t)+A/3*sin(3*f*t + phi).*exp(-a*t);
% y2=flipdim(-A*sin(f*t + phi).*exp(-a*t),2)+flipdim(-A/2*sin(2*f*t + phi).*exp(-a*t),2)+flipdim(-A/3*sin(3*f*t + phi).*exp(-a*t),2);
% y=[y2(1:end-1),y1];
% t=0:dt:2*dt*1e3;
% 
% % axis([0 1 -2.2 2.2]);
% 
% figure
% subplot(1,2,1);title('Time Trace','FontWeight','bold')
% plot(t,y)
% xlabel('Time (Sec)','FontWeight','bold')
% ylabel('Amplitude (a.u.)','FontWeight','bold')
% subplot(1,2,2);title('Power Spectrum','FontWeight','bold')
% freq=(0:length(t))/dt/length(t);
% ampl=(abs(fft(y))/length(t)).^2;
% plot(freq(1:end-1),ampl)
% xlabel('Frequency (Hz)','FontWeight','bold')
% ylabel('Amplitude (a.u.)','FontWeight','bold')

