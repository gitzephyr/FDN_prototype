%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Apr 10 23:07:41 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : Tue Apr 11 22:45:06 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- LPF and buffer -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% testing LPF and (circular) buffer delay line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LPFnBuff < audioPlugin
    properties 
        % LPF Coeff
        B0 = 0.985;
        % Coeff before and after del_buffer
        % cN = rand(1);
        % Input Gain
        Gain = 0.5;
        % LPF
        yLPFprev = [0 0];
        % init delay in samples
        Delay = 0.001;
    end
    properties (Access = private)
        % Delay Line
        z1 = zeros(220500,2); % 220500
        % index
        BufferIndex = 1;
        % Delay time in samples
        DSamples = 1;
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Gain','DisplayName','Dry','Mapping',{'lin',0,1}),...
            audioPluginParameter('B0','DisplayName','LPF Coeff','Mapping',{'lin',0,1}),...
            audioPluginParameter('Delay','DisplayName','Delay','Label','ms','Mapping',{'lin',0.001,1}));
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
                out(i,:) = (in(i,:) * plugin.Gain) + plugin.yLPFprev;
                % update delay line
                plugin.z1(writeIndex,:) = in(i,:);
                
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
        function set.Delay(plugin, val)
            plugin.Delay = val;
            plugin.DSamples = floor(getSampleRate(plugin)*val);
        end
    end
end
