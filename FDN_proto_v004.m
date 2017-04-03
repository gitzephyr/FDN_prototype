%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Fri Mar 19 14:30:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : Mon Apr  3 21:16:15 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% prototype 004 
% 16 TDL + 16 Delay + LOWPASS Filter
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
%              v  v   v   v   v
%              ----------------                           |
%              [      +       ]                           |
%              ----------------                           |
%                     |                                   |
%                     |                                   |
%                     -----------------------------------[+]----->y(n)
%
% [LATE REVERB] == FDN of 16 delays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pick a sound file
ls snd/
%% load a sound file
clear all; close all; clc;
[x, fs] = audioread('snd/singing.wav');
segm = zeros(fs,1);
x = [segm; x; segm];
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Opera voice');
%% or create an impulse
clear all; close all; clc;
fs = 44100;
x = zeros(1,5*fs);
x(1) = 1;
x = [zeros(1,fs), x, zeros(1,fs)];
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
plot(t,x); xlabel('Seconds'); ylabel('Amplitude');
title('Impulse');
%% or impulse of white noise
clear all; close all; clc;
fs = 44100;
len = fs*0.05;
wgNoise = wgn(len,1,0);
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
y = zeros(length(x),1);
b = 0.29*ones(1,16);
c = 0.4*ones(1,16);
% Gain coefficient |g|<1
g = 0.219999999999999999999;
%% Puckette Feedback Matrix 
% ??????????????????????
%% using Hadamard Matrix
A = g*(1/2)*hadamard(16);
%% Tapped Delay Line
M = zeros(1,16);
t = 0.004;
M(1) = fs*t;
for i = 2:16
    r = randi(40);
    M(i) = M(i-1) + r;
end
k = rand(1,16);
nX = zeros(size(x));
for n = sum(M)+1:length(x)
    nX(n) = k(1)*x(n) + k(2)*x(n-M(1)) + k(3)*x(n-M(2)) + k(4)*x(n-M(3)) + k(5)*x(n-M(4))...
        + k(6)*x(n-M(5)) + k(7)*x(n-M(6)) + k(8)*x(n-M(7)) + k(9)*x(n-M(8)) + k(10)*x(n-M(9))...
        + k(11)*x(n-M(10)) + k(12)*x(n-M(11)) + k(13)*x(n-M(12)) + k(14)*x(n-M(13))...
        + k(15)*x(n-M(14)) + k(16)*x(n-M(15));
end
%% 16 Delay lines, use prime
% m = [89 97 107 113 149 211 263 293 401 421 433 443 577 601 641 661]';
% m = [193 373 421 499 569 617 677 751 823 907 929 947 971 991 1019 1039]';
% m = [443 1949 4409 5417 6421 7537 8863 9049 10799 11177 12791 13679 14891 15287 16339 17657]';
m = DelayLineLengths(16,fs,0.05)
%% Delay lines
z1 = zeros(1,max(m));
z2 = zeros(1,max(m));
z3 = zeros(1,max(m));
z4 = zeros(1,max(m));
z5 = zeros(1,max(m));
z6 = zeros(1,max(m));
z7 = zeros(1,max(m));
z8 = zeros(1,max(m));
z9 = zeros(1,max(m));
z10 = zeros(1,max(m));
z11 = zeros(1,max(m));
z12 = zeros(1,max(m));
z13 = zeros(1,max(m));
z14 = zeros(1,max(m));
z15 = zeros(1,max(m));
z16 = zeros(1,max(m));
%% Loop
b0 = 0.3;
b1 = 1 - b0;
lastA = zeros(1,16);
for n = length(segm):length(y)
    temp = [z1(m(1)) z2(m(2)) z3(m(3)) z4(m(4)) z5(m(5)) z6(m(6)) z7(m(7)) z8(m(8))...
            z9(m(9)) z10(m(10)) z11(m(11)) z12(m(12)) z13(m(13)) z14(m(14)) z15(m(15)) z16(m(16))];
    
    y(n) = x(n) + c(1)*z1(m(1)) + c(2)*z2(m(2)) ...
        + c(3)*z3(m(3)) + c(4)*z4(m(4)) + c(5)*z5(m(5)) + c(6)*z6(m(6)) ...
        + c(7)*z7(m(7)) + c(8)*z8(m(8)) + c(9)*z9(m(9)) + c(10)*z10(m(10)) ...
        + c(11)*z11(m(11)) + c(12)*z12(m(12)) + c(13)*z13(m(13)) + c(14)*z14(m(14)) ...
        + c(15)*z15(m(15)) + c(16)*z16(m(16));
    
    % Lowpass filters
%     temp(1) = b1*temp(1) + b0*y(n-1);
%     temp(2) = b1*temp(2) + b0*y(n-1);
%     temp(3) = b1*temp(3) + b0*y(n-1);
%     temp(4) = b1*temp(4) + b0*y(n-1);
%     temp(5) = b1*temp(5) + b0*y(n-1);
%     temp(6) = b1*temp(6) + b0*y(n-1);
%     temp(7) = b1*temp(7) + b0*y(n-1);
%     temp(8) = b1*temp(8) + b0*y(n-1);
%     temp(9) = b1*temp(9) + b0*y(n-1);
%     temp(10) = b1*temp(10) + b0*y(n-1);
%     temp(11) = b1*temp(11) + b0*y(n-1);
%     temp(12) = b1*temp(12) + b0*y(n-1);
%     temp(13) = b1*temp(13) + b0*y(n-1);
%     temp(14) = b1*temp(14) + b0*y(n-1);
%     temp(15) = b1*temp(15) + b0*y(n-1);
%     temp(16) = b1*temp(16) + b0*y(n-1);
%     
%     z1 = [(x(n)*b(1) + temp*A(1,:)') z1(1:length(z1)-1)];
%     z2 = [(x(n)*b(2) + temp*A(2,:)') z2(1:length(z2)-1)]; 
%     z3 = [(x(n)*b(3) + temp*A(3,:)') z3(1:length(z3)-1)]; 
%     z4 = [(x(n)*b(4) + temp*A(4,:)') z4(1:length(z4)-1)];
%     z5 = [(x(n)*b(5) + temp*A(5,:)') z5(1:length(z5)-1)];
%     z6 = [(x(n)*b(6) + temp*A(6,:)') z6(1:length(z6)-1)]; 
%     z7 = [(x(n)*b(7) + temp*A(7,:)') z7(1:length(z7)-1)]; 
%     z8 = [(x(n)*b(8) + temp*A(8,:)') z8(1:length(z8)-1)];
%     z9 = [(x(n)*b(9) + temp*A(9,:)') z9(1:length(z9)-1)];
%     z10 = [(x(n)*b(10) + temp*A(10,:)') z10(1:length(z10)-1)]; 
%     z11 = [(x(n)*b(11) + temp*A(11,:)') z11(1:length(z11)-1)]; 
%     z12 = [(x(n)*b(12) + temp*A(12,:)') z12(1:length(z12)-1)];
%     z13 = [(x(n)*b(13) + temp*A(13,:)') z13(1:length(z13)-1)];
%     z14 = [(x(n)*b(14) + temp*A(14,:)') z14(1:length(z14)-1)]; 
%     z15 = [(x(n)*b(15) + temp*A(15,:)') z15(1:length(z15)-1)]; 
%     z16 = [(x(n)*b(16) + temp*A(16,:)') z16(1:length(z16)-1)];
    
    lastA(1) = b1*(temp*A(1,:)') + b0*lastA(1);
    lastA(2) = b1*(temp*A(2,:)') + b0*lastA(2);
    lastA(3) = b1*(temp*A(3,:)') + b0*lastA(3);
    lastA(4) = b1*(temp*A(4,:)') + b0*lastA(4);
    lastA(5) = b1*(temp*A(5,:)') + b0*lastA(5);
    lastA(6) = b1*(temp*A(6,:)') + b0*lastA(6);
    lastA(7) = b1*(temp*A(7,:)') + b0*lastA(7);
    lastA(8) = b1*(temp*A(8,:)') + b0*lastA(8);
    lastA(9) = b1*(temp*A(9,:)') + b0*lastA(9);
    lastA(10) = b1*(temp*A(10,:)') + b0*lastA(10);
    lastA(11) = b1*(temp*A(11,:)') + b0*lastA(11);
    lastA(12) = b1*(temp*A(12,:)') + b0*lastA(12);
    lastA(13) = b1*(temp*A(13,:)') + b0*lastA(13);
    lastA(14) = b1*(temp*A(14,:)') + b0*lastA(14);
    lastA(15) = b1*(temp*A(15,:)') + b0*lastA(15);
    lastA(16) = b1*(temp*A(16,:)') + b0*lastA(16);
    
    
    z1 = [(x(n)*b(1) + lastA(1)) z1(1:length(z1)-1)];
    z2 = [(x(n)*b(2) + lastA(2)) z2(1:length(z2)-1)]; 
    z3 = [(x(n)*b(3) + lastA(3)) z3(1:length(z3)-1)]; 
    z4 = [(x(n)*b(4) + lastA(4)) z4(1:length(z4)-1)];
    z5 = [(x(n)*b(5) + lastA(5)) z5(1:length(z5)-1)];
    z6 = [(x(n)*b(6) + lastA(6)) z6(1:length(z6)-1)]; 
    z7 = [(x(n)*b(7) + lastA(7)) z7(1:length(z7)-1)]; 
    z8 = [(x(n)*b(8) + lastA(8)) z8(1:length(z8)-1)];
    z9 = [(x(n)*b(9) + lastA(9)) z9(1:length(z9)-1)];
    z10 = [(x(n)*b(10) + lastA(10)) z10(1:length(z10)-1)]; 
    z11 = [(x(n)*b(11) + lastA(11)) z11(1:length(z11)-1)]; 
    z12 = [(x(n)*b(12) + lastA(12)) z12(1:length(z12)-1)];
    z13 = [(x(n)*b(13) + lastA(13)) z13(1:length(z13)-1)];
    z14 = [(x(n)*b(14) + lastA(14)) z14(1:length(z14)-1)]; 
    z15 = [(x(n)*b(15) + lastA(15)) z15(1:length(z15)-1)]; 
    z16 = [(x(n)*b(16) + lastA(16)) z16(1:length(z16)-1)];
end
%% Plot
dt = 1/fs;
t = 0:dt:(length(x)*dt)-dt;
figure(1)
subplot(4,1,1)
plot(t,x,'k'); 
subplot(4,1,2)
plot(t,nX,'g');
subplot(4,1,3)
plot(t,y)
subplot(4,1,4)
plot(t,nX); hold on;
plot(t,y);
xlabel('Seconds'); ylabel('Amplitude');
title('Feedback Delay Network');
% legend('reverb','original');
%% Spectrogram
figure(1) 
subplot(3,1,1);
title('original sound file');
specgram(x);
subplot(3,1,2);
specgram(nX);
subplot(3,1,3)
title('reverb');
specgram(y);
%% combine early reflections and late reflections
yy = y + nX;
%% 
soundsc(yy,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
