%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : March 15th, 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feedforward comb filter
% Echo simulator when b0 = 1 and bM = g (|bM| < 1)
% computational physical model of a single discrete echo.
% y(n) = b0*x(n) + bM*x(n-m)
clear all; close all; clc;
% [snd, sr] = audioread('anechoic_rec/b1-sc_44_1kHz.wav');
[snd, sr] = audioread('snd/singing.wav');
% soundsc(snd,sr);    % original sound
[s, c] = size(snd);
% N = samples ago
% N/sr = seconds ago
N = floor(sr*1.1);
del_snd = snd;  % set up a new array, same size as the old one
b0 = 1;
bM = 0.9;
for n = N+1:length(snd)
        del_snd(n) = b0*snd(n) + bM*snd(n-N);
end
plot(del_snd, 'k'); hold on;
plot(snd, 'g'); 
legend('filtered sound','original sound');
soundsc(del_snd,sr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feedback Comb Filter
% y(n) = b0*x(n) - aM*y(n-M)
clear all; close all; clc;
% [snd, sr] = audioread('anechoic_rec/b1-sc_44_1kHz.wav');
[snd, sr] = audioread('snd/singing.wav');
% soundsc(snd, sr);     % original sound
[s, c] = size(snd);
M = floor(sr*0.1);
del_snd = snd;          % set up a new array, same size as the old one
b0 = 0.5;
aM = 0.7;
for n = M+1:length(snd)
    del_snd(n) = b0*snd(n) - aM*del_snd(n-M);
end
plot(del_snd, 'k'); hold on;
plot(snd, 'g'); 
legend('filtered sound','original sound');
soundsc(del_snd,sr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computational physical model of a series of echoes!
% When y(n) = x(n) + g*y(n-M)
% Computational model of an ideal plane wave bouncing back and forth
% two parallel walls
clear all; close all; clc;
% [snd, sr] = audioread('anechoic_rec/b1-sc_44_1kHz.wav');
[snd, sr] = audioread('snd/singing.wav');
% soundsc(snd, sr);     % original sound
[s, c] = size(snd);
M = floor(sr*0.1);
del_snd = snd;          % set up a new array, same size as the old one
g = 0.6;
for n = M+1:length(snd)
    del_snd(n) = snd(n) + g*del_snd(n-M);
end
plot(del_snd, 'k'); hold on;
plot(snd, 'g'); 
legend('filtered sound','original sound');
soundsc(del_snd,sr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Version v02
% Sometimes y(n) is taken from the end of the delay line instead of the 
% beginning
% y(n) = bM*x(n-M) - aM*y(n-M)
clear all; close all; clc;
% [snd, sr] = audioread('anechoic_rec/b1-sc_44_1kHz.wav');
[snd, sr] = audioread('snd/singing.wav');
% soundsc(snd, sr);     % original sound
[s, c] = size(snd);
M = floor(sr*0.1);
del_snd = snd;          % set up a new array, same size as the old one
bM = 0.1;
aM = 0.5;
for n = M+1:length(snd)
    del_snd(n) = bM*snd(n-M) - aM*del_snd(n-M);
end
% plot(del_snd, 'k'); hold on;
% plot(snd, 'g'); 
% legend('filtered sound','original sound');
soundsc(del_snd,sr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMB FILTER
clear all; close all;
g1 = (0.5)^3;        % Some specific coefficients
g2 = (0.9)^5;
B = [1 0 0 g1];      % Feedforward coefficients, M1=3
A = [1 0 0 0 0 g2];  % Feedback coefficients, M2=5
N = 1000;            % Number of signal samples
x = rand(N,1);       % Random test input signal
y = filter(B,A,x);   % Matlab and Octave compatible
[y1,state] = filter(B,A,x);       % filter 1st block x1
[y2,state] = filter(B,A,x,state); % filter 2nd block x2
% The example coefficients,  g1 = 0.5^3 = 0.125$ and
% g2 = 0.9^5 = 0.59049, are chosen to place all filter zeros at radius 
% 0.5 and all filter poles at radius  0.9 in the complex z plane 
% (as we shall see below).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% efr.m - frequency response computation in Matlab/Octave
clear all; close all; clc;
% Example filter:
g1 = 0.5^3; B = [1 0 0 g1];      % Feedforward coeffs
g2 = 0.9^5; A = [1 0 0 0 0 g2];  % Feedback coefficients

Nfft = 1024;         % FFT size
Nspec = 1+Nfft/2;    % Show only positive frequencies
f=[0:Nspec-1]/Nfft;  % Frequency axis
Xnum = fft(B,Nfft);  % Frequency response of FIR part
Xden = fft(A,Nfft);  % Frequency response, feedback part
X = Xnum ./ Xden;    % Should check for divide by zero!
freqz(B,A,Nspec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%