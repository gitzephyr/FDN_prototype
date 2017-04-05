%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Mar 30 17:33:00 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : Wed Apr  5 20:13:37 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~ -*- Prime Power Delay-Line Lengths -*- ~~~~~~~~~~~~~~~~ %%
% DESCRIPTION
%   it is convenient to choose each delay-line length M_i, as an integer
%   power of a distinct prime number p_i. With this choice, the
%   delay-line lengths are always coprime.
%
%   Physical Audio Signal Processing
%   for Virtual Musical Instruments and Audio Effects
%   Julius O. Smith III
%   p. 113
%   
% INPUT
%   - N     : numbers of delay lines
%   - fs    : sampling frequency
%   - t     : first delay in ms
%
% OUTPUT
%   - m     : desired delay time in samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = DelayLineLengths(N,fs,t)
p = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131];
M = zeros(1,N);
M(1) = fs*t;
for i = 2:N
    r = randi(100);
    M(i) = M(i-1) + r;
end
% M = [fs*0.005 fs*0.009 fs*0.015 fs*0.020...
%      fs*0.028 fs*0.034 fs*0.040 fs*0.055 ...
%      fs*0.06 fs*0.065 fs*0.070 fs*0.077 ...
%      fs*0.08 fs*0.087 fs*0.099 fs*0.12];
M = round(M);
d = log(M(1:N));
n = log(p(1:N));
x = d./n;
nx = round(x);
m = p(1:N).^nx;
end

% del = p(1:4).^(round(log(M)./log(p(1:4))))
