%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Mon Apr 17 14:12:16 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : Thu Apr 20 22:57:43 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Head shadow model by Brown and Duda
% An Efficient Hrtf Model For 3-D Sound (1997)
% by C. Phillip Brown , Richard O. Duda
%% ~~~~~~~~~~~~~~~~~~~~~~ -*- Head Shadow -*- ~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
% Head shadow + Pinna Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef myHeadShadow < audioPlugin
    properties
        beta = (2*343)/0.0931;
        theta = 0;
        distance = 4;
    end
    properties (Access = private)
        thetaRad = 0;
        a = 0.0931;
        lastInL = 0;
        lastInR = 0;
        lastYL = 0;
        lastYR = 0;
        % Head Shadow buffers, indexes and samples
        bufferL = zeros(192001,1);
        bufferR = zeros(192001,1);
        bufferIndexL = 1;
        bufferIndexR = 1;
        NSamplesL = 1;
        NSamplesR = 1;
        
        % Pinna Model delays, buffers, indexes and samples
        % -- left
        bufferPM1L = zeros(192001,1);
        bufferPM2L = zeros(192001,1);
        bufferPM3L = zeros(192001,1);
        bufferPM4L = zeros(192001,1);
        bufferPM5L = zeros(192001,1);
        % -- rigth
        bufferPM1R = zeros(192001,1);
        bufferPM2R = zeros(192001,1);
        bufferPM3R = zeros(192001,1);
        bufferPM4R = zeros(192001,1);
        bufferPM5R = zeros(192001,1);
        % -- left
        bufferIndexPM1L = 1;
        bufferIndexPM2L = 1;
        bufferIndexPM3L = 1;
        bufferIndexPM4L = 1;
        bufferIndexPM5L = 1;
        % -- rigth
        bufferIndexPM1R = 1;
        bufferIndexPM2R = 1;
        bufferIndexPM3R = 1;
        bufferIndexPM4R = 1;
        bufferIndexPM5R = 1;
        % -- left
        NSamplesPM1L = 1;
        NSamplesPM2L = 1;
        NSamplesPM3L = 1;
        NSamplesPM4L = 1;
        NSamplesPM5L = 1;
        % -- right
        NSamplesPM1R = 1;
        NSamplesPM2R = 1;
        NSamplesPM3R = 1;
        NSamplesPM4R = 1;
        NSamplesPM5R = 1;
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('theta','DisplayName','theta','Label','Degree','Mapping',{'lin',-90,90}),...
            audioPluginParameter('distance','DisplayName','Distance','Mapping',{'lin',2,150}));
    end
    methods
        function out = process(plugin, in)
            out = zeros(size(in));
            fs = getSampleRate(plugin);
            t = 1/fs;
            tbeta = t*plugin.beta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- Head Shadow model
            % --- write indexes
            writeIndexL = plugin.bufferIndexL;
            writeIndexR = plugin.bufferIndexR;
            % --- read indexes
            readIndexL = writeIndexL - plugin.NSamplesL;
            readIndexR = writeIndexR - plugin.NSamplesR;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- Pinna Model
            % --- write indexes left and right
            wIndexPM1L = plugin.bufferIndexPM1L;
            wIndexPM2L = plugin.bufferIndexPM2L;
            wIndexPM3L = plugin.bufferIndexPM3L;
            wIndexPM4L = plugin.bufferIndexPM4L;
            wIndexPM5L = plugin.bufferIndexPM5L;
            
            wIndexPM1R = plugin.bufferIndexPM1R;
            wIndexPM2R = plugin.bufferIndexPM2R;
            wIndexPM3R = plugin.bufferIndexPM3R;
            wIndexPM4R = plugin.bufferIndexPM4R;
            wIndexPM5R = plugin.bufferIndexPM5R;
            % --- read indexes, left and right
            rIndexPM1L = wIndexPM1L - plugin.NSamplesPM1L;
            rIndexPM2L = wIndexPM2L - plugin.NSamplesPM2L;
            rIndexPM3L = wIndexPM3L - plugin.NSamplesPM3L;
            rIndexPM4L = wIndexPM4L - plugin.NSamplesPM4L;
            rIndexPM5L = wIndexPM5L - plugin.NSamplesPM5L;
            
            rIndexPM1R = wIndexPM1R - plugin.NSamplesPM1R;
            rIndexPM2R = wIndexPM2R - plugin.NSamplesPM2R;
            rIndexPM3R = wIndexPM3R - plugin.NSamplesPM3R;
            rIndexPM4R = wIndexPM4R - plugin.NSamplesPM4R;
            rIndexPM5R = wIndexPM5R - plugin.NSamplesPM5R;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- Head Shadow
            if readIndexL <= 0 
                readIndexL = readIndexL + 192001;
            end
            if readIndexR <= 0
                readIndexR = readIndexR + 192001;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- Pinna Model
            if rIndexPM1L <= 0
                rIndexPM1L = rIndexPM1L + 192001;
            end
            if rIndexPM2L <= 0
                rIndexPM2L = rIndexPM2L + 192001;
            end
            if rIndexPM3L <= 0
                rIndexPM3L = rIndexPM3L + 192001;
            end
            if rIndexPM4L <= 0
                rIndexPM4L = rIndexPM4L + 192001;
            end
            if rIndexPM5L <= 0
                rIndexPM5L = rIndexPM5L + 192001;
            end
            
            if rIndexPM1R <= 0
                rIndexPM1R = rIndexPM1R + 192001;
            end
            if rIndexPM2R <= 0
                rIndexPM2R = rIndexPM2R + 192001;
            end
            if rIndexPM3R <= 0
                rIndexPM3R = rIndexPM3R + 192001;
            end
            if rIndexPM4R <= 0
                rIndexPM4R = rIndexPM4R + 192001;
            end
            if rIndexPM5R <= 0
                rIndexPM5R = rIndexPM5R + 192001;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- process stuff
            for i = 1:size(in,1)
                % -- Head shadow
                plugin.bufferL(writeIndexL) = in(i,1); 
                plugin.bufferR(writeIndexR) = in(i,2);
                
                alpha_l = 1 - sin(plugin.theta*(pi/180));
                alpha_r = 1 + sin(plugin.theta*(pi/180)); 
                b0=2+tbeta; 
                b1=-2+tbeta;
                a0_l=2*alpha_l-tbeta; 
                a1_l=-2*alpha_l+tbeta;
                a0_r=2*alpha_r-tbeta;
                a1_r=-2*alpha_r+tbeta;
                % -- Head shadow 
                yl = ((a0_l*plugin.bufferL(readIndexL))+(a1_l*plugin.lastInL))-(b1*plugin.lastYL)*(1/b0);
                yr = ((a0_r*plugin.bufferR(readIndexR))+(a1_r*plugin.lastInR))-(b1*plugin.lastYR)*(1/b0);
                % -- Pinna Model
                plugin.bufferPM1L(wIndexPM1L) = yl;
                plugin.bufferPM2L(wIndexPM2L) = yl;
                plugin.bufferPM3L(wIndexPM3L) = yl;
                plugin.bufferPM4L(wIndexPM4L) = yl;
                plugin.bufferPM5L(wIndexPM5L) = yl;
                
                plugin.bufferPM1R(wIndexPM1R) = yr;
                plugin.bufferPM2R(wIndexPM2R) = yr;
                plugin.bufferPM3R(wIndexPM3R) = yr;
                plugin.bufferPM4R(wIndexPM4R) = yr;
                plugin.bufferPM5R(wIndexPM5R) = yr;
                %
                xPL = 0.5*plugin.bufferPM1L(rIndexPM1L) + ...
                      -1*plugin.bufferPM2L(rIndexPM2L) + ...
                      0.5*plugin.bufferPM3L(rIndexPM3L) + ...
                      -0.25*plugin.bufferPM4L(rIndexPM4L) + ...
                      0.25*plugin.bufferPM5L(rIndexPM5L);
                  
                xPR = 0.5*plugin.bufferPM1R(rIndexPM1R) + ...
                      -1*plugin.bufferPM2R(rIndexPM2R) + ...
                      0.5*plugin.bufferPM3R(rIndexPM3R) + ...
                      -0.25*plugin.bufferPM4R(rIndexPM4R) + ...
                      0.25*plugin.bufferPM5R(rIndexPM5R);
                
                out(i,1) = (yl + xPL)*(2/plugin.distance);
                out(i,2) = (yr + xPR)*(2/plugin.distance);
                
                plugin.lastInL = plugin.bufferL(readIndexL);
                plugin.lastInR = plugin.bufferR(readIndexR);
                plugin.lastYL = yl;
                plugin.lastYR = yr;
                
                % -- updating indexes
                % --- head shadow
                writeIndexL = writeIndexL + 1;
                writeIndexR = writeIndexR + 1;
                % --- Pinna Model
                wIndexPM1L = wIndexPM1L + 1;
                wIndexPM2L = wIndexPM2L + 1;
                wIndexPM3L = wIndexPM3L + 1;
                wIndexPM4L = wIndexPM4L + 1;
                wIndexPM5L = wIndexPM5L + 1;
                
                wIndexPM1R = wIndexPM1R + 1;
                wIndexPM2R = wIndexPM2R + 1;
                wIndexPM3R = wIndexPM3R + 1;
                wIndexPM4R = wIndexPM4R + 1;
                wIndexPM5R = wIndexPM5R + 1;
                
                if writeIndexL > 192001
                    writeIndexL = 1;
                end
                if writeIndexR > 192001
                    writeIndexR = 1;
                end
                % -- left
                if wIndexPM1L > 192001
                    wIndexPM1L = 1;
                end
                if wIndexPM2L > 192001
                    wIndexPM2L = 1;
                end
                if wIndexPM3L > 192001
                    wIndexPM3L = 1;
                end
                if wIndexPM4L > 192001
                    wIndexPM4L = 1;
                end
                if wIndexPM5L > 192001
                    wIndexPM5L = 1;
                end
                % -- right
                if wIndexPM1R > 192001
                    wIndexPM1R = 1;
                end
                if wIndexPM2R > 192001
                    wIndexPM2R = 1;
                end
                if wIndexPM3R > 192001
                    wIndexPM3R = 1;
                end
                if wIndexPM4R > 192001
                    wIndexPM4R = 1;
                end
                if wIndexPM5R > 192001
                    wIndexPM5R = 1;
                end
                % -- head shadow
                readIndexL = readIndexL + 1;
                readIndexR = readIndexR + 1;
                
                if readIndexL > 192001
                    readIndexL = 1;
                end
                if readIndexR > 192001
                    readIndexR = 1;
                end
                % -- Pinna model
                rIndexPM1L = rIndexPM1L + 1;
                rIndexPM2L = rIndexPM2L + 1;
                rIndexPM3L = rIndexPM3L + 1;
                rIndexPM4L = rIndexPM4L + 1;
                rIndexPM5L = rIndexPM5L + 1;
                
                rIndexPM1R = rIndexPM1R + 1;
                rIndexPM2R = rIndexPM2R + 1;
                rIndexPM3R = rIndexPM3R + 1;
                rIndexPM4R = rIndexPM4R + 1;
                rIndexPM5R = rIndexPM5R + 1;
                
                if rIndexPM1L > 192001
                    rIndexPM1L = 1;
                end
                if rIndexPM2L > 192001
                    rIndexPM2L = 1;
                end
                if rIndexPM3L > 192001
                    rIndexPM3L = 1;
                end
                if rIndexPM4L > 192001
                    rIndexPM4L = 1;
                end
                if rIndexPM5L > 192001
                    rIndexPM5L = 1;
                end
                
                if rIndexPM1R > 192001
                    rIndexPM1R = 1;
                end
                if rIndexPM2R > 192001
                    rIndexPM2R = 1;
                end
                if rIndexPM3R > 192001
                    rIndexPM3R = 1;
                end
                if rIndexPM4R > 192001
                    rIndexPM4R = 1;
                end
                if rIndexPM5R > 192001
                    rIndexPM5R = 1;
                end
                
            end
            plugin.bufferIndexL = writeIndexL;
            plugin.bufferIndexR = writeIndexR;
            
            plugin.bufferIndexPM1L = wIndexPM1L;
            plugin.bufferIndexPM2L = wIndexPM2L;
            plugin.bufferIndexPM3L = wIndexPM3L;
            plugin.bufferIndexPM4L = wIndexPM4L;
            plugin.bufferIndexPM5L = wIndexPM5L;
            
            plugin.bufferIndexPM1R = wIndexPM1R;
            plugin.bufferIndexPM2R = wIndexPM2R;
            plugin.bufferIndexPM3R = wIndexPM3R;
            plugin.bufferIndexPM4R = wIndexPM4R;
            plugin.bufferIndexPM5R = wIndexPM5R;
        end
        function set.theta(plugin, val)
            plugin.theta = val;
            plugin.thetaRad = val*(pi/180);
            % -- Pinna Echoes
            plugin.NSamplesPM1L = floor(1*cos((plugin.thetaRad*-1)/2)*sin(1*(1.57-0))+2);
            plugin.NSamplesPM2L = floor(5*cos((plugin.thetaRad*-1)/2)*sin(0.5*(1.57-0))+4);
            plugin.NSamplesPM3L = floor(5*cos((plugin.thetaRad*-1)/2)*sin(0.5*(1.57-0))+7);
            plugin.NSamplesPM4L = floor(5*cos((plugin.thetaRad*-1)/2)*sin(0.5*(1.57-0))+11);
            plugin.NSamplesPM5L = floor(5*cos((plugin.thetaRad*-1)/2)*sin(0.5*(1.57-0))+13);

            plugin.NSamplesPM1R = floor(1*cos(plugin.thetaRad/2)*sin(0.5*(1.57-0))+2);
            plugin.NSamplesPM2R = floor(5*cos(plugin.thetaRad/2)*sin(0.5*(1.57-0))+4);
            plugin.NSamplesPM3R = floor(5*cos(plugin.thetaRad/2)*sin(0.5*(1.57-0))+5);
            plugin.NSamplesPM4R = floor(5*cos(plugin.thetaRad/2)*sin(0.5*(1.57-0))+7);
            plugin.NSamplesPM5R = floor(5*cos(plugin.thetaRad/2)*sin(0.5*(1.57-0))+13);
            % -- Head Shadow 
            plugin.NSamplesL = floor((plugin.a-plugin.a*sin(plugin.thetaRad)/343)*getSampleRate(plugin));
            plugin.NSamplesR = floor((plugin.a+plugin.a*(plugin.thetaRad)/343)*getSampleRate(plugin));
            
%             ptr = 'L-------------------------------------------------------'
%             plugin.NSamplesPM1L
%             plugin.NSamplesPM2L
%             plugin.NSamplesPM3L
%             plugin.NSamplesPM4L
%             plugin.NSamplesPM5L
%             ptr = 'R-------------------------------------------------------'
%             plugin.NSamplesPM1R
%             plugin.NSamplesPM2R
%             plugin.NSamplesPM3R
%             plugin.NSamplesPM4R
%             plugin.NSamplesPM5R
        end
    end
end
