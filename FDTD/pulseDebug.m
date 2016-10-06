dt = 1E-6;
t = 0:dt:0.2;
f0 = 60E3;
BW = 0.3;
offset = 0;
plotPulse = 1;

x = pulse( t, dt, f0, BW, offset, 0, plotPulse);