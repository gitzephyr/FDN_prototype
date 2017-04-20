%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Thu Apr 20 16:32:55 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% prototype 005 
% 4 TDL going through an HRTF model + 4 Delay + LOWPASS Filter
% based on: 
% - Physical Audio Signal Processing
%   for Virtual Musical Instruments and Audio Effects
%   Julius O. Smith III
% p. 65-67, p.85-127
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structure:
% x(n)------->[TAPPED DELAY LINE]------>[LATE REVERB]------
%              |  |   |   |   |                           |
%              |  |   |   |   |                           |
%              v  v   v   v   v                           |
%              ----------------                           |
%              [     HRTF     ]                           |
%              ----------------                           |
%              |  |   |   |   |                           |
%              |  |   |   |   |                           |
%              v  v   v   v   v                           |
%              ----------------                           |
%              [      +       ]                           |
%              ----------------                           |
%                     |                                   |
%                     |                                   |
%                     -----------------------------------[+]----->y(n)
%
% [LATE REVERB] >> FDN of 4 delays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pick a sound file
ls snd/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load a sound file
clear all; close all; clc;
[x, fs] = audioread('snd/singing.wav');
segm = zeros(fs,1);
x = [segm; x; segm];
dt = 1/fs;
tt = 0:dt:(length(x)*dt)-dt;
plot(tt,x); xlabel('Seconds'); ylabel('Amplitude');
title('Opera voice');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% play it!
soundsc(x,fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize parameters
nfdn = 4;
y = zeros(length(x),1);
% b = 0.29*ones(1,16);
b = rand(1,nfdn);
% c = 0.4*ones(1,16);
c = rand(1,nfdn);
% Gain coefficient |g|<1
g = 0.219999999999999999999;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using Hadamard Matrix
A = g*(1/2)*hadamard(nfdn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tapped Delay Line
M = zeros(1,nfdn);
t = 0.009;
M(1) = fs*t;
for i = 2:16
    r = randi(100);
    M(i) = M(i-1) + r;
end
k = rand(1,nfdn);
nX = zeros(size(x));
for n = sum(M)+1:length(x)
    nX(n) = k(1)*x(n) + k(2)*x(n-M(1)) + k(3)*x(n-M(2)) + k(4)*x(n-M(3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% play early reflection 
soundsc(nX,fs)
%% HRTF
yl = zeros(1,length(x));
yr = zeros(1,length(x));
std_circ = 58.5; % male, cm. From wikipedia
% std_circ = 53; % female, cm
a = std_circ/(2*pi); % cm
a = a/100; % meters
beta = (2*343)/a; 
t = 1/fs;
tbeta = (t*beta);
n = 1;
dT = 45;
for theta = dT
    alpha_l=1-sin(theta*(pi/180));
    alpha_r=1+sin(theta*(pi/180));

    % filter coefficients
    b0=2+tbeta; 
    b1=-2+tbeta;
    a0_l=2*alpha_l-tbeta;
    a1_l=-2*alpha_l+tbeta;
    a0_r=2*alpha_r-tbeta;
    a1_r=-2*alpha_r+tbeta;
    
    for i = 2:length(nX)
        yl(i)=((a0_l*nX(i))+(a1_l*nX(i-1)))-(b1*yl(i-1))*(1/b0);
        yr(i)=((a0_r*nX(i))+(a1_r*nX(i-1)))-(b1*yr(i-1))*(1/b0);
    end
end
%% merge yLeft and yRigth
yHrtf = [yl; yr];
%% plot it!
plot(yl,fs,'m'); hold on;
plot(yr,fs,'b')
%% play it!
soundsc(yHrtf,fs)
%% 4 Delay lines, use prime
m = prime_power_delays(fs,nfdn,10,20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDN Delay lines
z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main FDN Loop
b0 = 0.3;
b1 = 1 - b0;
lastA = zeros(1,16);
for n = length(segm):length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4))];
    
    y(n) = x(n) + c(1)*z1(m(1)) + c(2)*z2(m(2)) + c(3)*z3(m(3)) + c(4)*z4(m(4));
    
    lastA(1) = b1*(temp*A(1,:)') + b0*lastA(1);
    lastA(2) = b1*(temp*A(2,:)') + b0*lastA(2);
    lastA(3) = b1*(temp*A(3,:)') + b0*lastA(3);
    lastA(4) = b1*(temp*A(4,:)') + b0*lastA(4);
    
    z1 = [(x(n)*b(1) + lastA(1)) z1(1:length(z1)-1)];
    z2 = [(x(n)*b(2) + lastA(2)) z2(1:length(z2)-1)]; 
    z3 = [(x(n)*b(3) + lastA(3)) z3(1:length(z3)-1)]; 
    z4 = [(x(n)*b(4) + lastA(4)) z4(1:length(z4)-1)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% combine early reflections and late reflections
yy = y + nX;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
figure(1)
subplot(5,1,1)
plot(t,x,'k'); 
xlabel('sec'); ylabel('amp');
title('INPUT');
subplot(5,1,2)
plot(t,nX,'g');
xlabel('sec'); ylabel('amp');
title('Early reflection');
subplot(5,1,3)
plot(t,y)
xlabel('sec'); ylabel('amp');
title('Late reflections');
subplot(5,1,4)
plot(t,nX); hold on;
plot(t,y);
xlabel('sec'); ylabel('amp');
title('Early + Late Reflection');
subplot(5,1,5)
plot(t,yy)
xlabel('sec'); ylabel('amp');
title('OUTPUT');
% legend('reverb','original');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spectrogram
figure(1) 
subplot(4,1,1);
specgram(x);
subplot(4,1,2);
specgram(nX);
subplot(4,1,3)
specgram(y);
subplot(4,1,4)
specgram(yy);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
soundsc(yy,fs)
%% 
soundsc(y,fs)
%% 
soundsc(x,fs)
%%
soundsc(nX,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


