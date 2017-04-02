%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Fri Mar 19 14:30:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : 
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
        B0 = 0.97;
        % Delay Time
        Delay = 0.5;
        % Coeff before and after del_buffer
        B = 0.3;
        C = 0.6;
        % Feedback Matrix
        Dampening = 0.4;
        % Input Gain
        Gain = 0.5;
        % LPF
        yLast = 0;
    end
    properties (Access = private)
        % Delay Lines
        z1 = zeros(220500,2); % 220500
        z2 = zeros(220500,2);
        z3 = zeros(220500,2);
        z4 = zeros(220500,2);
        % 
        BufferIndex = 1;
        % Delay times
        NSamples = [32 243 625 343]';
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...  %<---
            audioPluginParameter('Gain',...         %<---
            'DisplayName','Input Gain',...           %<---
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening',...
            'DisplayName','Dampening',...
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('B',...
            'DisplayName','B Coeff',...
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('C',...
            'DisplayName','C Coeff',...
            'Mapping',{'lin',0,1}));
%         PluginInterface = audioPluginInterface(...
%             audioPluginParameter('Delay',...
%             'DisplayName','Delay',...
%             'Label','seconds',...
%             'Mapping',{'lin',0,5}))
%             audioPluginParameter('B0',...
%             'DisplayName','LPF Coeff',...
%             'Mapping',{'lin',0,1})

    end
    methods
        function out = process(plugin, in)
            % b and c coff - buffers
            bN = plugin.B*ones(1,4);
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
%                 temp(1) = B1*temp(1) + plugin.B0*plugin.yLast;
%                 temp(2) = B1*temp(2) + plugin.B0*plugin.yLast;
%                 temp(3) = B1*temp(3) + plugin.B0*plugin.yLast;
%                 temp(4) = B1*temp(4) + plugin.B0*plugin.yLast;
%                 
%                 plugin.yLast = y(1);
                
                plugin.z1(writeIndex,:) = in(i,:)*bN(1) + temp*A(1,:)';
                plugin.z2(writeIndex,:) = in(i,:)*bN(2) + temp*A(2,:)';
                plugin.z3(writeIndex,:) = in(i,:)*bN(3) + temp*A(3,:)';
                plugin.z4(writeIndex,:) = in(i,:)*bN(4) + temp*A(4,:)';
                
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
%         function set.Delay(plugin, val)
%             plugin.Delay = val;
%             plugin.NSamples = floor(getSampleRate(plugin)*val);
%         end
    end
end