%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Julius O. Smith (jos at ccrma.stanford.edu)
% Created on        : Wed Apr 5 18:13:26 CEST 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : Wed Apr  5 19:12:37 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- prime power delays -*- ~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% Prime Power Delay Line Lengths
%
% USAGE:
%   prime_power_delays(fs,N,pathmin,pathmax)
%
% WHERE
%   N = positive integer up to 16
%       (for higher powers of 2, extend 'primes' array below.)
%   pathmin = minimum acoustic ray length in the reverberator (in meters)
%   pathmax = maximum acoustic ray length (meters) - think "room size"
%
% REFERENCE:
%   https://ccrma.stanford.edu/~jos/pasp/Prime_Power_Delay_Line.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = prime_power_delays(fs,N,pathmin,pathmax)
    Np = N;
    i = [1:Np];
    pathmin = 55;
    pathmax = 63;
    prime = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131];

    % Prime Power Bounds [matlab: floor(log(maxdel)./log(primes(53)))]
    maxdel=8192; % more than 63 meters at 44100 samples/sec & 343 m/s
    % ppbs = [13,8,5,4,3,3,3,3,2,2,2,2,2,2,2,2]; % 8192 is enough for all
    % ppb(i) = take(i+1,ppbs);

    % Approximate desired delay-line lengths using powers of distinct primes:
    c = 343; % soundspeed in m/s at 20 degrees C for dry air
    dmin = fs*pathmin/c;
    dmax = fs*pathmax/c;
    dl = dmin * (dmax/dmin).^(i/(Np-1)); % desired delay in samples
    ppwr = floor(log(dl)./log(prime(1:Np))); % best prime power
    m = prime(1:Np).^ppwr; % each delay a power of a distinct prime
end
