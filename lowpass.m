%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Mar 20 23:12:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~~~~~~ -*- LOWPASS FILTER -*- ~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% lowpass.m
%% LOWPASS FILTER - Simple - white noise
clear all; close all;
fs = 44100;
sec = 2;
len = fs*sec;
a = -0.25; b = 0.25;
% interval [a, b]
segm = zeros(fs,1);
wnoise = a + (b-a).*rand(len,1);
wnoise = [segm; wnoise; segm];

figure(1)
subplot(2,1,1);
plot(wnoise);
xlabel('Seconds'); ylabel('Amplitude');
title('Low Pass Filter');

y = zeros(len,1);
for i = fs:length(wnoise)
    y(i) = wnoise(i) + wnoise(i-1);
end
subplot(2,1,2);
plot(y); xlabel('Seconds'); ylabel('Amplitude');
title('Low Pass Filter');
%% WARNING!! if you're using headphones, turn down the volume
soundsc(wnoise,fs)
%% WARNING!! if you're using headphones, turn down the volume
soundsc(y,fs)
%% Spectrogram
figure(1) 
subplot(2,1,1);
specgram(wnoise);
subplot(2,1,2);
specgram(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOWPASS FILTER - Simple - sound file
clear all; close all;
% [x, fs] = audioread('../anechoic_rec/singing.wav');
[x, fs] = audioread('snd/singing.wav');
segm = zeros(fs,1);
x = [segm; x; segm];
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;

len = length(x);
y = zeros(len,1);
g = 0.1;
for i = fs:length(x)
    % y(i) = g*x(i) + g*x(i-1) + g*x(i-2) + g*x(i-3) + g*x(i-4) + g*x(i-5) + g*x(i-5) + g*x(i-6);
    y(i) = g*x(i) + (1-g)*y(i-1);
    % y(i) = x(i) + x(i-1) ;
end

figure(1)
subplot(2,1,1);
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Singing - Opera voice');
subplot(2,1,2);
plot(y); xlabel('Seconds'); ylabel('Amplitude');
title('Low Pass Filter');
soundsc(x,fs);
pause(round(length(x)/fs))
soundsc(y,fs);
%% Spectrogram
figure(2)
subplot(2,1,1);
specgram(x);
subplot(2,1,2);
specgram(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATLAB IMPLEMENTATION WITH FILTER FUNCTION
clear all; close all; clc;
%% White noise
fs = 44100;
a = -1; b = 1;
% interval [a, b]
len = fs*2;
x = a + (b-a).*rand(len,1);
%% sound file
[x, fs] = audioread('snd/singing.wav');
%%
segm = zeros(fs,1);
x = [segm; x; segm];
B = [0.5 0.5]; 
figure(1);
freqz(B,1);
%%
figure(2);
zplane(B,1);
y = filter(B,1,x);
%% Plot waveform
figure(1)
subplot(2,1,1);
plot(x); xlabel('Seconds'); ylabel('Amplitude');
title('Singing - Opera voice');
subplot(2,1,2);
plot(y); xlabel('Seconds'); ylabel('Amplitude');
title('Low Pass Filter');
%% Spectrogram
figure(2) 
subplot(2,1,1);
specgram(x);
subplot(2,1,2);
specgram(y);
%%
soundsc(x,fs);
%%
soundsc(y,fs)