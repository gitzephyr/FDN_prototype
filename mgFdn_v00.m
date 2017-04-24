%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Thu Apr 20 22:01:19 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : Mon Apr 24 22:27:33 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~~~~~~~~~~~ -*- mgFdn -*- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% structure:
% x(n)------->[TAPPED DELAY LINE]------>[LATE REVERB]------
%              |  |   |   |   |                           |
%              |  |   |   |   |                           |
%              v  v   v   v   v                           |
%              ----------------                     ---------------
%              [ HEAD SHADOW  ]                     [ LATE REVERB ]
%              ----------------                     [    FDN      ]
%              |  |   |   |   |                     ---------------
%              |  |   |   |   |                           |
%              |  |   |   |   |                           |
%              v  v   v   v   v                           |
%              ----------------                           |
%              [  PINNA MODEL ]                           |
%              ----------------                           |
%               |           |                             |
%               |           |                             |
%               |           v                             v
%               v           -----------------------------[+]----->y_left(n)
%               -----------------------------------------[+]----->y_right(n)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N Tapped Delay Line
% Head Shadow > on each TDL, two channels as output.
% Pinna Model > 
% [LATE REVERB] >> FDN of 4 delays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
classdef mgFdn_v00 < audioPlugin
    properties
        % LPF coeff
        B0 = 0.1;
        % Fdn B and C gain
        rB = rand(1,4);
        rC = rand(1,4);
        B = 0.5;
        C = 0.5;
        % Feedback Matrix
        Dampening = 0.1;
        hadamardMtrx = hadamard(4);
        % Gain
        Gain = 0.5;
        % init pathmin and pathmax
        pathmin = 3;
        pathmax = 5;
        % temp variables
        dmin = 0;
        dmax = 0;
        % sound speed
        c = 343;
        % prime numbers needed for delay lines
        prime = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53];
        % TDL coeff
        k = rand(1,4);
        kEarly = 0.5;
        % 
        beta = (2*343)/0.0931;
        theta = 0;
        distance = 4;
    end
    properties (Access = private)
        % TDL
        tdl = zeros(192001,2);
        % FDN Delay Lines
        z1 = zeros(192001,2);
        z2 = zeros(192001,2);
        z3 = zeros(192001,2);
        z4 = zeros(192001,2);
        % Head Shadow
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
        % Index
        BufferIndex = 1;
        BufferTapIndex = 1;
        % Tapped Delay Lines 
        % TSamples = [343 432 464 476 635 683 747 802 867 922 995 1048 1148 1170 1181 1192];
        TSamples = [1549 1613 1699 1764]; 
        % FDN Delay times
        NSamples = [512 729 625 343 1331 2197 289 361 529 841 961 1369 1681 1849 2209 2809];
        % NSamples = [1681 1849 2209 2809];
        % Last output
        lpfPrev = zeros(1,4);
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('theta','DisplayName','theta','Label','Degree','Mapping',{'lin',-90,90}),...
            audioPluginParameter('distance','DisplayName','Distance','Mapping',{'lin',2,150}),...
            audioPluginParameter('Gain','DisplayName','Dry','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('C','DisplayName','Wet','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('kEarly','DisplayName','Pre-delay','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening','DisplayName','Dampening','Label','%','Mapping',{'lin',0,0.5}),...
            audioPluginParameter('B0','DisplayName','LowPass Filter','Mapping',{'lin',0,1}),...
            audioPluginParameter('pathmin','DisplayName','Min Room Size','Label','meters','Mapping',{'int',1,50}),...
            audioPluginParameter('pathmax','DisplayName','Max Room Size','Label','meters','Mapping',{'int',1,50}));
    end
    methods
        function out = process(plugin, in)
            out = zeros(size(in));
            fs = getSampleRate(plugin);
            t = 1/fs;
            tbeta = t*plugin.beta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- write indexes 
            wIndex = plugin.BufferIndex;
            wIndexTap = plugin.BufferTapIndex;
            wIndexL = plugin.bufferIndexL;
            wIndexR = plugin.bufferIndexR;
            % -- read indexes
            rIndexZ1 = wIndex - plugin.NSamples(1);
            rIndexZ2 = wIndex - plugin.NSamples(2);
            rIndexZ3 = wIndex - plugin.NSamples(3);
            rIndexZ4 = wIndex - plugin.NSamples(4);
            rIndexTDL_1 = wIndexTap - plugin.TSamples(1);
            rIndexTDL_2 = wIndexTap - plugin.TSamples(2);
            rIndexTDL_3 = wIndexTap - plugin.TSamples(3);
            rIndexTDL_4 = wIndexTap - plugin.TSamples(4);
            rIndexL = wIndexL - plugin.NSamplesL;
            rIndexR = wIndexR - plugin.NSamplesR;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % wrap it
            if rIndexZ1 <= 0
                rIndexZ1 = rIndexZ1 + 192001;
            end
            if rIndexZ2 <= 0
                rIndexZ2 = rIndexZ2 + 192001;
            end
            if rIndexZ3 <= 0
                rIndexZ3 = rIndexZ3 + 192001;
            end
            if rIndexZ4 <= 0
                rIndexZ4 = rIndexZ4 + 192001;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rIndexTDL_1 <= 0
                rIndexTDL_1 = rIndexTDL_1 + 192001;
            end
            if rIndexTDL_2 <= 0
                rIndexTDL_2 = rIndexTDL_2 + 192001;
            end
            if rIndexTDL_3 <= 0
                rIndexTDL_3 = rIndexTDL_3 + 192001;
            end
            if rIndexTDL_4 <= 0
                rIndexTDL_4 = rIndexTDL_4 + 192001;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rIndexL <= 0 
                rIndexL = rIndexL + 192001;
            end
            if rIndexR <= 0
                rIndexR = rIndexR + 192001;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -- processing
            for i = 1:size(in,1)
                
                
                % B and C coeff FDN
                BB = plugin.rB*plugin.B;
                CC = plugin.rC*plugin.C;
                % k coeff TDL
                KK = plugin.k*plugin.kEarly;
                % LPF coeff
                B1 = 1 - plugin.B0;
                % Feedback Matrix
                A = plugin.Dampening*(1/2)*plugin.hadamardMtrx;
                % Head shadow 
                alpha_l = 1 - sin(plugin.theta*(pi/180));
                alpha_r = 1 + sin(plugin.theta*(pi/180)); 
                b0=2+tbeta; 
                b1=-2+tbeta;
                a0_l=2*alpha_l-tbeta; 
                a1_l=-2*alpha_l+tbeta;
                a0_r=2*alpha_r-tbeta;
                a1_r=-2*alpha_r+tbeta;
                % ---------------------------------------------------------
                % tapped delay line
                tdl_out = KK(1)*plugin.tdl(rIndexTDL_1,:) + ...
                          KK(2)*plugin.tdl(rIndexTDL_2,:) + ...
                          KK(3)*plugin.tdl(rIndexTDL_3,:) + ...
                          KK(4)*plugin.tdl(rIndexTDL_4,:);
                % ---------------------------------------------------------
                % -- Head shadow 
                yl = ((a0_l*plugin.bufferL(rIndexL))+(a1_l*plugin.lastInL))-(b1*plugin.lastYL)*(1/b0);
                yr = ((a0_r*plugin.bufferR(rIndexR))+(a1_r*plugin.lastInR))-(b1*plugin.lastYR)*(1/b0);
                
                plugin.lastInL = plugin.bufferL(rIndexL);
                plugin.lastInR = plugin.bufferR(rIndexR);
                plugin.lastYL = yl;
                plugin.lastYR = yr;
                % ---------------------------------------------------------
                % temp delay samples
                temp = [plugin.z1(rIndexZ1) plugin.z2(rIndexZ2)...
                        plugin.z3(rIndexZ3) plugin.z4(rIndexZ4)];    
                % ---------------------------------------------------------
                % LowPass Filters after each delay line
                plugin.lpfPrev(1) = plugin.B0*plugin.z1(rIndexZ1) + B1*plugin.lpfPrev(1);
                plugin.lpfPrev(2) = plugin.B0*plugin.z2(rIndexZ2) + B1*plugin.lpfPrev(2);
                plugin.lpfPrev(3) = plugin.B0*plugin.z3(rIndexZ3) + B1*plugin.lpfPrev(3);
                plugin.lpfPrev(4) = plugin.B0*plugin.z4(rIndexZ4) + B1*plugin.lpfPrev(4);
                % ---------------------------------------------------------
                % equation -- output
                out(i,1) = yl + CC(1)*plugin.lpfPrev(1) + ...
                            CC(2)*plugin.lpfPrev(2) + CC(3)*plugin.lpfPrev(3) + ...
                            CC(4)*plugin.lpfPrev(4);
                        
                out(i,2) = yr + CC(1)*plugin.lpfPrev(1) + ...
                            CC(2)*plugin.lpfPrev(2) + CC(3)*plugin.lpfPrev(3) + ...
                            CC(4)*plugin.lpfPrev(4);
                % ---------------------------------------------------------
                % -- write to circular buffer
                plugin.tdl(wIndexTap,:) = in(i,:) * plugin.Gain;
                plugin.bufferL(wIndexL) = tdl_out(1,1);
                plugin.bufferR(wIndexR) = tdl_out(1,2);
                
                plugin.z1(wIndex,:) = plugin.tdl(rIndexTDL_1,:)*BB(1) + temp*A(1,:)';
                plugin.z2(wIndex,:) = plugin.tdl(rIndexTDL_2,:)*BB(2) + temp*A(2,:)';
                plugin.z3(wIndex,:) = plugin.tdl(rIndexTDL_3,:)*BB(3) + temp*A(3,:)';
                plugin.z4(wIndex,:) = plugin.tdl(rIndexTDL_4,:)*BB(4) + temp*A(4,:)';
                % ---------------------------------------------------------
                % update indexes
                wIndex = wIndex + 1;
                wIndexTap = wIndexTap + 1;
                wIndexL = wIndexL + 1;
                wIndexR = wIndexR + 1;
                
                rIndexZ1 = rIndexZ1 + 1;
                rIndexZ2 = rIndexZ2 + 1;
                rIndexZ3 = rIndexZ3 + 1;
                rIndexZ4 = rIndexZ4 + 1;
                rIndexTDL_1 = rIndexTDL_1 + 1;
                rIndexTDL_2 = rIndexTDL_2 + 1;
                rIndexTDL_3 = rIndexTDL_3 + 1;
                rIndexTDL_4 = rIndexTDL_4 + 1;
                rIndexL = rIndexL + 1;
                rIndexR = rIndexR + 1;
                % ---------------------------------------------------------
                if wIndex > 192001
                    wIndex = 1;
                end
                if wIndexTap > 192001
                    wIndexTap = 1;
                end
                if wIndexL > 192001
                    wIndexL = 1;
                end
                if wIndexR > 192001
                    wIndexR = 1;
                end
                % ---------------------------------------------------------
                if rIndexZ1 > 192001
                    rIndexZ1 = 1;
                end
                if rIndexZ2 > 192001
                    rIndexZ2 = 1;
                end
                if rIndexZ3 > 192001
                    rIndexZ3 = 1;
                end
                if rIndexZ4 > 192001
                    rIndexZ4 = 1;
                end
                if rIndexTDL_1 > 192001
                    rIndexTDL_1 = 1;
                end
                if rIndexTDL_2 > 192001
                    rIndexTDL_2 = 1;
                end
                if rIndexTDL_3 > 192001
                    rIndexTDL_3 = 1;
                end
                if rIndexTDL_4 > 192001
                    rIndexTDL_4 = 1;
                end
                if rIndexL > 192001
                    rIndexL = 1;
                end
                if rIndexR > 192001
                    rIndexR = 1;
                end
            end
            plugin.BufferIndex = wIndex;
            plugin.BufferTapIndex = wIndexTap;
            plugin.bufferIndexL = wIndexL;
            plugin.bufferIndexR = wIndexR;
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
            ppwr = floor(0.5 + log(dl)./log(plugin.prime(i))); % best prime power
            plugin.NSamples(i) = plugin.prime(i).^ppwr; % each delay a power of a distinct prime
%             ptr = '--------------------------------------------------------'
%             plugin.NSamples(i)
%             plugin.TSamples(i)
%             plugin.NSamplesL
%             plugin.NSamplesR
        end
        %------------------------------------------------------------------
        function set.pathmax(plugin, val)
            Np = 4;
            i = [1:Np];
            plugin.pathmax = val;
            plugin.dmax = getSampleRate(plugin)*val/plugin.c;
            dl = plugin.dmin * (plugin.dmax/plugin.dmin).^(i/(Np-1));
            ppwr = floor(0.5 + log(dl)./log(plugin.prime(i)));
            plugin.NSamples(i) = plugin.prime(i).^ppwr;
%             ptr = '--------------------------------------------------------'
%             plugin.NSamples(i)
%             plugin.TSamples(i) 
%             plugin.NSamplesL
%             plugin.NSamplesR
        end
        function set.theta(plugin, val)
            plugin.theta = val;
            plugin.thetaRad = val*(pi/180);
            % -- Head Shadow 
            plugin.NSamplesL = floor((plugin.a-plugin.a*sin(plugin.thetaRad)/343)*getSampleRate(plugin))-4000;
            plugin.NSamplesR = floor((plugin.a+plugin.a*(plugin.thetaRad)/343)*getSampleRate(plugin))-4000;
%             ptr = '--------------------------------------------------------'
%             plugin.NSamplesL
%             plugin.NSamplesR
        end
    end
end
