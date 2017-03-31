%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : T. Lokki
% Created on        : Fri Mar 17 12:19:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% prototype 000
%% pick a sound file
ls snd/
%%
clear all; close all; clc;
[x, fs] = audioread('snd/singing.wav');
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Opera voice');
%% play it!
soundsc(x,fs);
%% Init parameters
y = zeros(1,length(x));
b = 0.4*ones(1,8);
c = 0.4*ones(1,8);
% Gain coefficient |g|<1
g = 0.969999999;
%% Puckette Feedback Matrix 
a = [0 1 1 0;
    -1 0 0 -1;
    1 0 0 -1;
    0 1 -1 0];
A = g*(1/sqrt(2))*a
%% Hadamard Matrix
A = g*(1/2)*hadamard(4);
%% Delay, use prime
% m =[149 211 263 293]';
% m = [577 601 641 661]';
m = DelayLineLengths(4)
%% Delay lines
z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
%% Loop FDN
for n = 1:length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4))]; 
    y(n) = x(n) + c(1)*z1(m(1)) + c(2)*z2(m(2)) ...
        + c(3)*z3(m(3)) + c(4)*z4(m(4));
    z1 = [(x(n)*b(1) + temp*A(1,:)') z1(1:length(z1)-1)];
    z2 = [(x(n)*b(2) + temp*A(2,:)') z2(1:length(z2)-1)]; 
    z3 = [(x(n)*b(3) + temp*A(3,:)') z3(1:length(z3)-1)]; 
    z4 = [(x(n)*b(4) + temp*A(4,:)') z4(1:length(z4)-1)]; 
end
%% Plot
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); hold on;
plot(t,y); xlabel('Seconds'); ylabel('Amplitude');
title('Feedback Delay Network');
%% Spectrogram
figure(1) 
subplot(2,1,1);
specgram(x);
subplot(2,1,2);
specgram(y);
%%
soundsc(y,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

