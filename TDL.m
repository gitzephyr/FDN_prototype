%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Mar 15 13:22:00 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tapped Delay Line (TDL)
% 
clear all; close all; clc;
[x, sr] = audioread('snd/singing.wav');
% soundsc(x, sr);                   % original sound
[s, c] = size(x);
M3 = floor(sr*1.2);                   % total delay line length
nTDL = 4;                           % number of delay lines
step = M3/nTDL;
M = zeros(1,nTDL);                  % vector of delays in samples
M(nTDL) = M3;
for n = 1:nTDL-1
    M(n) = floor(step*n);
end
b = [0.3 0.2 0.2 0.2 0.9];              % vector of b0, bM1..bMn
offset = zeros(sum(M), 1);             % generate a zeros segment concat to x
x = [offset; x; offset];
y = zeros(size(x));                 % set up a new array of zeros!, same size as x
for n = sum(M)+1:length(x)
    y(n) = b(1)*x(n) + b(2)*x(n-M(1)) + b(3)*x(n-M(2)) + b(4)*x(n-M(3)) + b(5)*x(n-M(4));
end
%% 
subplot(2,1,1)
plot(y, 'k');
legend('filtered sound');
subplot(2,1,2)
plot(x, 'g'); 
legend('original sound');
soundsc(y, sr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tapped Delay Line (TDL) v002
% smaller delay lines
clear all; close all; clc;
[x, sr] = audioread('snd/singing.wav');
% soundsc(x, sr);                   % original sound
[s, c] = size(x);
M = [sr*0.002 sr*0.009 sr*0.030 sr*0.1];
b = [0.3 0.2 0.2 0.2 0.4];              % vector of b0, bM1..bMn
offset = zeros(sum(M), 1);             % generate a zeros segment concat to x
x = [offset; x; offset];
y = zeros(size(x));                 % set up a new array of zeros!, same size as x
for n = sum(M)+1:length(x)
    y(n) = b(1)*x(n) + b(2)*x(n-M(1)) + b(3)*x(n-M(2)) + b(4)*x(n-M(3)) + b(5)*x(n-M(4));
end
%% 
subplot(2,1,1)
plot(y, 'k');
legend('filtered sound');
subplot(2,1,2)
plot(x, 'g'); 
legend('original sound');
soundsc(y, sr);