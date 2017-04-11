%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Fri Mar 19 14:30:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : Mon Apr  3 21:17:25 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
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
classdef myFDN < audioPlugin
    properties 
 % LPF Coeff
        B0 = 0.985;
        % Delay Time
        % Delay = 1;
        % Coeff before and after del_buffer
        bFdn = rand(1,4);
        cFdn = rand(1);
        % Feedback Matrix
        Dampening = 0.1;
        % Input Gain
        Gain = 0.5;
        % LPF
        yLast = 0;
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
        % tapped delay line coeff k, random
        k = rand(1,16);
        % wet 
        Wet = 0.5;
    end
    properties (Access = private)
        % Delay Lines
        z1 = zeros(220500,2);
        z2 = zeros(220500,2);
        z3 = zeros(220500,2);
        z4 = zeros(220500,2);
        % FDN index
        BufferIndex = 1;
        % Delay times
        NSamples = zeros(1,4);
        % Last output of LPF
        lpfPrev = zeros(1,16);
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Gain','DisplayName','Dry','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening','DisplayName','Dampening','Mapping',{'lin',0,0.5}),...
            audioPluginParameter('cFdn','DisplayName','C Coeff','Mapping',{'lin',0,1}),...
            audioPluginParameter('B0','DisplayName','LPF Coeff','Mapping',{'lin',0,1}),...
            audioPluginParameter('pathmin','DisplayName','RoomSizeMin','Label','meters','Mapping',{'int',1,50}),...
            audioPluginParameter('pathmax','DisplayName','RoomSizeMax','Label','meters','Mapping',{'int',1,50}));
    end
    methods
        function out = process(plugin, in)
            % b and c coff - buffers
            cN = plugin.cFdn*ones(1,4);
            
            B1 = 1 - plugin.B0;
            
            out = zeros(size(in));
            writeIndex = plugin.BufferIndex;
            
            readIndexZ1 = writeIndex - plugin.NSamples(1);
            readIndexZ2 = writeIndex - plugin.NSamples(2);
            readIndexZ3 = writeIndex - plugin.NSamples(3);
            readIndexZ4 = writeIndex - plugin.NSamples(4);
            
            if readIndexZ1 <= 0
                readIndexZ1 = readIndexZ1 + 220500;
            end
            if readIndexZ2 <= 0
                readIndexZ2 = readIndexZ2 + 220500;
            end
            if readIndexZ3 <= 0
                readIndexZ3 = readIndexZ3 + 220500;
            end
            if readIndexZ4 <= 0
                readIndexZ4 = readIndexZ4 + 220500;
            end
            
            for i = 1:size(in,1)
                % lpf coeff
                B1 = 1 - plugin.B0;
                % feedback matrix, controls diffusion
                A = plugin.Dampening*(1/2)*hadamard(4);
                % tmp delay samples
                temp = [plugin.z1(readIndexZ1) plugin.z2(readIndexZ2)...
                    plugin.z3(readIndexZ3) plugin.z4(readIndexZ4)];
                % LowPass Filters after each delay line
                plugin.lpfPrev(1) = plugin.B0*plugin.z1(readIndexZ1) + B1*plugin.lpfPrev(1);
                plugin.lpfPrev(2) = plugin.B0*plugin.z2(readIndexZ2) + B1*plugin.lpfPrev(2);
                plugin.lpfPrev(3) = plugin.B0*plugin.z3(readIndexZ3) + B1*plugin.lpfPrev(3);
                plugin.lpfPrev(4) = plugin.B0*plugin.z4(readIndexZ4) + B1*plugin.lpfPrev(4);
                % equation
%                 outBuffers = cN(1)*plugin.z1(Z1_readIndex) + ...
%                     cN(2)*plugin.z2(Z2_readIndex) + cN(3)*plugin.z3(Z3_readIndex) + ...
%                     cN(4)*plugin.z4(Z4_readIndex);
                outBuffers = cN(1) * plugin.lpfPrev(1) + cN(2) * plugin.lpfPrev(2) + ...
                    cN(3) * plugin.lpfPrev(3) + cN(4) * plugin.lpfPrev(4);
                y = (in(i,:) * plugin.Gain) + outBuffers;
                % output
                out(i,:) = y;
                
                % LOWPASS Filter
%                 plugin.lastA(1) = B1*(temp*A(1,:)') + plugin.B0*plugin.lastA(1);
%                 plugin.lastA(2) = B1*(temp*A(2,:)') + plugin.B0*plugin.lastA(2);
%                 plugin.lastA(3) = B1*(temp*A(3,:)') + plugin.B0*plugin.lastA(3);
%                 plugin.lastA(4) = B1*(temp*A(4,:)') + plugin.B0*plugin.lastA(4);
%                 temp(1) = B1*temp(1) + plugin.B0*plugin.yLast;
%                 temp(2) = B1*temp(2) + plugin.B0*plugin.yLast;
%                 temp(3) = B1*temp(3) + plugin.B0*plugin.yLast;
%                 temp(4) = B1*temp(4) + plugin.B0*plugin.yLast;
%                 plugin.yLast = sum(y)/2;
                
%                 plugin.z1(writeIndex,:) = plugin.lastA(1)*plugin.bN(1) + temp*A(1,:)';
%                 plugin.z2(writeIndex,:) = plugin.lastA(2)*plugin.bN(2) + temp*A(2,:)';
%                 plugin.z3(writeIndex,:) = plugin.lastA(3)*plugin.bN(3) + temp*A(3,:)';
%                 plugin.z4(writeIndex,:) = plugin.lastA(4)*plugin.bN(4) + temp*A(4,:)';
                % y
                plugin.z1(writeIndex,:) = in(i,:)*plugin.bFdn(1) + temp*A(1,:)';
                plugin.z2(writeIndex,:) = in(i,:)*plugin.bFdn(2) + temp*A(2,:)';
                plugin.z3(writeIndex,:) = in(i,:)*plugin.bFdn(3) + temp*A(3,:)';
                plugin.z4(writeIndex,:) = in(i,:)*plugin.bFdn(4) + temp*A(4,:)';
                
                writeIndex = writeIndex + 1;
                
                if writeIndex > 220500
                    writeIndex = 1;
                end
                
                readIndexZ1 = readIndexZ1 + 1;
                readIndexZ2 = readIndexZ2 + 1;
                readIndexZ3 = readIndexZ3 + 1;
                readIndexZ4 = readIndexZ4 + 1;
                
                if readIndexZ1 > 220500
                    readIndexZ1 = 1;
                end
                if readIndexZ2 > 220500
                    readIndexZ2 = 1;
                end
                if readIndexZ3 > 220500
                    readIndexZ3 = 1;
                end
                if readIndexZ4 > 220500
                    readIndexZ4 = 1;
                end
            end
            plugin.BufferIndex = writeIndex;
        end
        %------------------------------------------------------------------
        function set.pathmin(plugin, val)
            fs = getSampleRate(plugin);
            Np = 4;
            i = [1:Np];
            % Approximate desired delay-line lengths using powers of distinct primes:
            % c = 343; % soundspeed in m/s at 20 degrees C for dry air
            plugin.pathmin = val;
            plugin.dmin = fs*val/plugin.c;
            dl = plugin.dmin * (plugin.dmax/plugin.dmin).^(i/(Np-1)); % desired delay in samples
            ppwr = floor(log(dl)./log(plugin.prime(1:Np))); % best prime power
            plugin.NSamples = plugin.prime(1:Np).^ppwr; % each delay a power of a distinct prime
        end
        %------------------------------------------------------------------
        function set.pathmax(plugin, val)
            Np = 4;
            i = [1:Np];
            plugin.pathmax = val;
            plugin.dmax = getSampleRate(plugin)*val/plugin.c;
            dl = plugin.dmin * (plugin.dmax/plugin.dmin).^(i/(Np-1));
            ppwr = floor(0.5 + log(dl)./log(plugin.prime(1:Np)));
            plugin.NSamples = plugin.prime(1:Np).^ppwr; 
        end
    end
end
