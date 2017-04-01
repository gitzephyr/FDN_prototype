%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Mar 20 09:57:00 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% FDN_proto_v002.m
% 8 Delay + LOWPASS Filter
%% pick a sound file
ls snd/
%% load a sound file
clear all; close all; clc;
[x, fs] = audioread('snd/singing.wav');
segm = zeros(fs,1);
x = [segm; x; segm; segm];
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Sound file');
%% or create an impulse
clear all; close all; clc;
fs = 44100;
x = zeros(1,5*fs);
x(1) = 1;
segm = zeros(1,fs);
x = [segm, x, segm];
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Impulse');
%% or impulse of white noise
clear all; close all; clc;
fs = 44100;
len = fs*0.1;
wgNoise = wgn(4410,1,0);
m  =  ((max(wgNoise) - 1) / 10);
x = wgNoise * m;
if (max(x) >= 1)
    x = x * 0.89;
end
x = x.*hamming(len);
segm = zeros(fs,1);
x = [segm; x; segm; segm; segm];
plot(x); xlabel('Seconds'); ylabel('Amplitude');
title('Impulse');
%% play it!
soundsc(x,fs);
%%
y = zeros(1,length(x));
b = 0.5*ones(1,8);
c = 0.5*ones(1,8);
% Gain coefficient |g|<1
g = 0.269999999999999999999;
%% Puckette and Stautner Matrix
% ??????
%% using Hadamard Matrix
A = g*(1/2)*hadamard(8);
%% 8 Delay lines, use prime
% m = [89 97 107 113 149 211 263 293]';
% m = [401 421 433 443 577 601 641 661]' ;
% m = [89 263 443 661 829 1039 1423 1693]';
% m = [443 829 1039 1423 1693 2467 3371 4409]';
% m = [443 761 1321 1949 3121 4969 5827 7537]';
% m = [443 1949 4409 9049 11177 12791 15287 17657]';
% m = [4409 11177 17137 26347 29629 35117 37619 40013]';
m = DelayLineLengths(8);
%% Delay lines
z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
z5 = zeros(1,max(m));
z6 = zeros(1,max(m));
z7 = zeros(1,max(m));
z8 = zeros(1,max(m));
%% Loop
b0 = 0.98;
b1 = 1 - b0;
for n = length(segm):length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4)) z5(m(5)) z6(m(6)) z7(m(7)) z8(m(8))];
    
    y(n) = x(n) + c(1)*z1(m(1)) + c(2)*z2(m(2)) ...
        + c(3)*z3(m(3)) + c(4)*z4(m(4)) + c(5)*z5(m(5)) + c(6)*z6(m(6)) ...
        + c(7)*z7(m(7)) + c(8)*z8(m(8));
    
    % Lowpass filters
    temp(1) = b1*temp(1) + b0*y(n-1);
    temp(2) = b1*temp(2) + b0*y(n-1);
    temp(3) = b1*temp(3) + b0*y(n-1);
    temp(4) = b1*temp(4) + b0*y(n-1);
    temp(5) = b1*temp(5) + b0*y(n-1);
    temp(6) = b1*temp(6) + b0*y(n-1);
    temp(7) = b1*temp(7) + b0*y(n-1);
    temp(8) = b1*temp(8) + b0*y(n-1);
    
    z1 = [(x(n)*b(1) + temp*A(1,:)') z1(1:length(z1)-1)];
    z2 = [(x(n)*b(2) + temp*A(2,:)') z2(1:length(z2)-1)]; 
    z3 = [(x(n)*b(3) + temp*A(3,:)') z3(1:length(z3)-1)]; 
    z4 = [(x(n)*b(4) + temp*A(4,:)') z4(1:length(z4)-1)];
    z5 = [(x(n)*b(5) + temp*A(5,:)') z5(1:length(z5)-1)];
    z6 = [(x(n)*b(6) + temp*A(6,:)') z6(1:length(z6)-1)]; 
    z7 = [(x(n)*b(7) + temp*A(7,:)') z7(1:length(z7)-1)]; 
    z8 = [(x(n)*b(8) + temp*A(8,:)') z8(1:length(z8)-1)];
end
%% Plot
dt = 1/fs; 
t = 0:dt:(length(x)*dt)-dt;
plot(t,y,'g'); hold on;
plot(t,x,'k'); xlabel('Seconds'); ylabel('Amplitude');
title('Feedback Delay Network');
legend('reverb','original');
%% Spectrogram
figure(1) 
subplot(2,1,1);
specgram(x);
subplot(2,1,2);
specgram(y);
%%
soundsc(y,fs),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%