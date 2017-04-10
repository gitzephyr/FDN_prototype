%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Apr 10 23:07:41 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- LPF and buffer -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% Real-time implementation of FDN
% you need the Audio System Toolbox!! 
% 
% 4 Delay + LOWPASS Filter
% based on: 
% - Physical Audio Signal Processing
%   for Virtual Musical Instruments and Audio Effects
%   Julius O. Smith III
% p. 65-67, p.85-127
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LPFnBuff < audioPlugin
    properties 
        % LPF Coeff
        B0 = 0.985;
        % Coeff before and after del_buffer
        cN = rand(1);
        % Input Gain
        Gain = 0.5;
        % LPF
        yLPFprev = [0 0];
        % init pathmin and pathmax
        pathmin = 3;
        pathmax = 5;
        % temp variables
        dmin = 0;
        dmax = 0;
        % sound speed
        c = 343;
        % prime numbers needed for delay lines
        prime = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131];
    end
    properties (Access = private)
        % Delay Line
        z1 = zeros(220500,2); % 220500
        % index
        BufferIndex = 1;
        % Delay time in samples
        DSamples = 44;
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Gain','DisplayName','Dry','Mapping',{'lin',0,1}),...
            audioPluginParameter('B0','DisplayName','LPF Coeff','Mapping',{'lin',0,1}),...
            audioPluginParameter('pathmin','DisplayName','RoomSizeMin','Label','meters','Mapping',{'int',1,50}),...
            audioPluginParameter('pathmax','DisplayName','RoomSizeMax','Label','meters','Mapping',{'int',1,50}));
    end
    methods
        function out = process(plugin, in)
            % init output
            out = zeros(size(in));
            % init index
            writeIndex = plugin.BufferIndex;
            readIndexZ = writeIndex - plugin.DSamples;
            % check if readIndex exceed bounds
            if readIndexZ <= 0
                readIndexZ = readIndexZ + 220500;
            end
            % 
            for i = 1:size(in,1)                
                B1 = 1 - plugin.B0;
                % LowPass Filters after delay line
                plugin.yLPFprev = plugin.B0*plugin.z1(readIndexZ) + B1*plugin.yLPFprev;
                % output - equation                
                out(i,:) = (in(i,:) * plugin.Gain) + plugin.cN * plugin.yLPFprev;
                % update delay line
                plugin.z1(writeIndex,:) = plugin.yLPFprev;
                
                % update writeIndex
                writeIndex = writeIndex + 1;
                % wrap it around
                if writeIndex > 220500
                    writeIndex = 1;
                end
                % update readIndex
                readIndexZ = readIndexZ + 1;
                % wrap it around
                if readIndexZ > 220500
                    readIndexZ = 1;
                end
            end
            plugin.BufferIndex = writeIndex;
        end
        %------------------------------------------------------------------
        function set.pathmin(plugin, val)
            Np = 1;
            i = [1:Np];
            % Approximate desired delay-line lengths using powers of distinct primes:
            % c = 343; % soundspeed in m/s at 20 degrees C for dry air
            plugin.pathmin = val;
            plugin.dmin = getSampleRate(plugin)*val/plugin.c;
            dl = plugin.dmin * (plugin.dmax/plugin.dmin).^(i/(Np-1)); % desired delay in samples
            ppwr = floor(log(dl)./log(plugin.prime(1:Np))); % best prime power
            plugin.DSamples = plugin.prime(1:Np).^ppwr; % each delay a power of a distinct prime
        end
        %------------------------------------------------------------------
        function set.pathmax(plugin, val)
            Np = 1;
            i = [1:Np];
            plugin.pathmax = val;
            plugin.dmax = getSampleRate(plugin)*val/plugin.c;
            dl = plugin.dmin * (plugin.dmax/plugin.dmin).^(i/(Np-1));
            ppwr = floor(0.5 + log(dl)./log(plugin.prime(1:Np)));
            plugin.DSamples = plugin.prime(1:Np).^ppwr; 
        end
    end
end
