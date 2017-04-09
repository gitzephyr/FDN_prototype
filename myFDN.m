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
        bN = rand(1,4);
        C = rand(1);
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
        z1 = zeros(220500,2); % 220500
        z2 = zeros(220500,2);
        z3 = zeros(220500,2);
        z4 = zeros(220500,2);
        % FDN index
        BufferIndex = 1;
        % Delay times
        NSamples = zeros(1,4);
        % Last output
        lastA = zeros(1,16);
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Gain','DisplayName','Dry','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening','DisplayName','Dampening','Mapping',{'lin',0,0.5}),...
            audioPluginParameter('C','DisplayName','C Coeff','Mapping',{'lin',0,1}),...
            audioPluginParameter('B0','DisplayName','LPF Coeff','Mapping',{'lin',0,1}),...
            audioPluginParameter('pathmin','DisplayName','RoomSizeMin','Label','meters','Mapping',{'int',1,50}),...
            audioPluginParameter('pathmax','DisplayName','RoomSizeMax','Label','meters','Mapping',{'int',1,50}));
    end
    methods
        function out = process(plugin, in)
            % b and c coff - buffers
            cN = plugin.C*ones(1,4);
            
%             B1 = 1 - plugin.B0;
            
            out = zeros(size(in));
            writeIndex = plugin.BufferIndex;
            
            Z1_readIndex = writeIndex - plugin.NSamples(1);
            Z2_readIndex = writeIndex - plugin.NSamples(2);
            Z3_readIndex = writeIndex - plugin.NSamples(3);
            Z4_readIndex = writeIndex - plugin.NSamples(4);
            
            if Z1_readIndex <= 0
                Z1_readIndex = Z1_readIndex + 220500;
            end
            if Z2_readIndex <= 0
                Z2_readIndex = Z2_readIndex + 220500;
            end
            if Z3_readIndex <= 0
                Z3_readIndex = Z3_readIndex + 220500;
            end
            if Z4_readIndex <= 0
                Z4_readIndex = Z4_readIndex + 220500;
            end
            
            for i = 1:size(in,1)
                
                B1 = 1 - plugin.B0;
                % feedback matrix
                % controls the diffusion
                % Feedback == Dampening
                A = plugin.Dampening*(1/2)*hadamard(4);
                
                temp = [plugin.z1(Z1_readIndex) plugin.z2(Z2_readIndex)...
                    plugin.z3(Z3_readIndex) plugin.z4(Z4_readIndex)];
                
                % equation
                y = (in(i,:) * plugin.Gain) + cN(1)*plugin.z1(Z1_readIndex) + ...
                    cN(2)*plugin.z2(Z2_readIndex) + cN(3)*plugin.z3(Z3_readIndex) + ...
                    cN(4)*plugin.z4(Z4_readIndex);
                % out
                out(i,:) = y;
                
                % LOWPASS Filter
%                 plugin.lastA(1) = B1*(temp*A(1,:)') + plugin.B0*plugin.lastA(1);
%                 plugin.lastA(2) = B1*(temp*A(2,:)') + plugin.B0*plugin.lastA(2);
%                 plugin.lastA(3) = B1*(temp*A(3,:)') + plugin.B0*plugin.lastA(3);
%                 plugin.lastA(4) = B1*(temp*A(4,:)') + plugin.B0*plugin.lastA(4);
                temp(1) = B1*temp(1) + plugin.B0*plugin.yLast;
                temp(2) = B1*temp(2) + plugin.B0*plugin.yLast;
                temp(3) = B1*temp(3) + plugin.B0*plugin.yLast;
                temp(4) = B1*temp(4) + plugin.B0*plugin.yLast;
                plugin.yLast = sum(y)/2;
                
                plugin.z1(writeIndex,:) = in(i,:)*plugin.bN(1) + temp*A(1,:)';
                plugin.z2(writeIndex,:) = in(i,:)*plugin.bN(2) + temp*A(2,:)';
                plugin.z3(writeIndex,:) = in(i,:)*plugin.bN(3) + temp*A(3,:)';
                plugin.z4(writeIndex,:) = in(i,:)*plugin.bN(4) + temp*A(4,:)';
                
%                 plugin.z1(writeIndex,:) = in(i,:)*bN(1) + plugin.lastA(1);
%                 plugin.z2(writeIndex,:) = in(i,:)*bN(2) + plugin.lastA(2);
%                 plugin.z3(writeIndex,:) = in(i,:)*bN(3) + plugin.lastA(3);
%                 plugin.z4(writeIndex,:) = in(i,:)*bN(4) + plugin.lastA(4);
                
                writeIndex = writeIndex + 1;
                
                if writeIndex > 220500
                    writeIndex = 1;
                end
                
                Z1_readIndex = Z1_readIndex + 1;
                Z2_readIndex = Z2_readIndex + 1;
                Z3_readIndex = Z3_readIndex + 1;
                Z4_readIndex = Z4_readIndex + 1;
                
                if Z1_readIndex > 220500
                    Z1_readIndex = 1;
                end
                if Z2_readIndex > 220500
                    Z2_readIndex = 1;
                end
                if Z3_readIndex > 220500
                    Z3_readIndex = 1;
                end
                if Z4_readIndex > 220500
                    Z4_readIndex = 1;
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
