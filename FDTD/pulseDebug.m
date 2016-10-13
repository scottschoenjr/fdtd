clear all
close all
clc


dt = 5E-7;
t = 0:dt:0.2;
f0 = 100E3;
BW = 0.5;
offset = 20E-3;
plotPulse = 1;

x = pulse( t, f0, BW, offset, 0, plotPulse);