%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Fri Mar 19 14:30:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% prototype 001 
% 4 Delay + LOWPASS Filter
%% pick a sound file
ls snd/
%%
clear all; close all; clc;
[x, fs] = audioread('snd/singing.wav');
segm = zeros(fs,1);
x = [segm; x; segm];
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Opera voice');
%% play it!
soundsc(x,fs);
%%
y = zeros(1,length(x));
b = 0.3*ones(1,8);
c = 0.6*ones(1,8);
% Gain coefficient |g|<1
g = 0.4;
%% Puckette Feedback Matrix 
a = [0 1 1 0;
    -1 0 0 -1;
    1 0 0 -1;
    0 1 -1 0];
A = (g*(1/sqrt(2)))*a;
%% using Hadamard Matrix
A = g*(1/2)*hadamard(4);
%% 4 Delay lines, use prime
% m =[149 211 263 293]';
% m = [401 421 433 443]';
% long
% m = [577 601 641 661]';
% short
% m = [89 97 107 113]';
% very long
% m = [919 941 971 997];
% tmp
% m = [89 443 1423 4409]';
% m = [1024 2187 3125 16807]';
m = DelayLineLengths(4);
%% Delay lines
z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
%% Loop
b0 = 0.25;
b1 = 1 - b0;
for n = length(segm):length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4))]; 
    
    y(n) = x(n) + c(1)*z1(m(1)) + c(2)*z2(m(2)) ...
        + c(3)*z3(m(3)) + c(4)*z4(m(4));
    
    % Lowpass filters
    temp(1) = b1*temp(1) + b0*y(n-1);
    temp(2) = b1*temp(2) + b0*y(n-1);
    temp(3) = b1*temp(3) + b0*y(n-1);
    temp(4) = b1*temp(4) + b0*y(n-1);
    
    % update buffers
    z1 = [(x(n)*b(1) + temp*A(1,:)') z1(1:length(z1)-1)];
    z2 = [(x(n)*b(2) + temp*A(2,:)') z2(1:length(z2)-1)]; 
    z3 = [(x(n)*b(3) + temp*A(3,:)') z3(1:length(z3)-1)]; 
    z4 = [(x(n)*b(4) + temp*A(4,:)') z4(1:length(z4)-1)]; 
end
%% Plot
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,y,'k'); hold on;
plot(t,x,'g'); xlabel('Seconds'); ylabel('Amplitude');
title('Feedback Delay Network');
legend('reverb','original');
%% Spectrogram
figure(1) 
subplot(2,1,1);
title('original sound file');
specgram(x);
subplot(2,1,2);
title('reverb');
specgram(y);
%%
soundsc(y,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

