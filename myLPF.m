%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Sun Apr 9 22:48:08 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~~~~~~~~ -*- LowPass Filter -*- ~~~~~~~~~~~~~~~~~~~~~~~ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef myLPF < audioPlugin
    properties 
        g = 0.1;
    end
    properties (Access = private)
        yLast = [0 0];
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('g','DisplayName','LPF Coeff','Mapping',{'lin',0,1}));
    end
    methods
        function out = process(plugin, in)
            out = zeros(size(in));
            for i = 1:size(in,1)
                tmp = plugin.g*in(i,:) + (1-plugin.g)*plugin.yLast;
                out(i,:) = tmp;
                plugin.yLast = tmp;
            end
        end
    end
end
