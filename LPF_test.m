clear all; close all; clc;
% [snd, sr] = audioread('anechoic_rec/b1-sc_44_1kHz.wav');
[snd, sr] = audioread('snd/singing.wav');
a = 0.1;
% y(n) = ax(n) + (1-a)x(n-1)
% feedforward 
B = [a (1-a)];
% feedback
A = [1 0];
freqz(B,A)
%%
a = 0.9;
% y(n) = ax(n) + (1-a)y(n-1)
% feedforward 
B = [a 1];
% feedback
A = [1 (1-a)];
freqz(B,A)
%%
y = filter(B,A,snd);
soundsc(y,sr)