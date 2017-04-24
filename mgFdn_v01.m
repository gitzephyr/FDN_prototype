%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Sun Apr 23 17:54:58 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : Mon Apr 24 22:27:37 CEST 2017
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
% 6 Tapped Delay Line
% Head Shadow > on each TDL, two channels as output.
% Pinna Model > 
% [LATE REVERB] >> FDN of 16 delays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
classdef mgFdn_v01 < audioPlugin
    properties
        % LPF coeff
        B0 = 0.1;
        % Fdn B and C gain
        rB = rand(1,16);
        rC = rand(1,16);
        B = 0.5;
        C = 0.5;
        % Feedback Matrix
        Dampening = 0.1;
        hadamardMtrx = hadamard(16);
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
        k = rand(1,6);
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
        z5 = zeros(192001,2);
        z6 = zeros(192001,2);
        z7 = zeros(192001,2);
        z8 = zeros(192001,2);
        z9 = zeros(192001,2);
        z10 = zeros(192001,2);
        z11 = zeros(192001,2);
        z12 = zeros(192001,2);
        z13 = zeros(192001,2);
        z14 = zeros(192001,2);
        z15 = zeros(192001,2);
        z16 = zeros(192001,2);        
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
        TSamples = [995 1181 1549 1613 1699 1764]; 
        % FDN Delay times
        NSamples = [512 729 625 343 1331 2197 289 361 529 841 961 1369 1681 1849 2209 2809];
        % Last output
        lpfPrev = zeros(1,16);
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('theta','DisplayName','theta','Label','Degree','Mapping',{'lin',-90,90}),...
            audioPluginParameter('Gain','DisplayName','Dry','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('C','DisplayName','Wet','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('distance','DisplayName','Distance','Mapping',{'lin',2,150}),...
            audioPluginParameter('kEarly','DisplayName','Pre-delay','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening','DisplayName','Dampening','Label','%','Mapping',{'lin',0,0.5}),...
            audioPluginParameter('B0','DisplayName','LowPass Filter','Mapping',{'lin',0,1}),...
            audioPluginParameter('pathmax','DisplayName','Room Size','Label','meters','Mapping',{'int',1,50}));
        % audioPluginParameter('pathmin','DisplayName','Min Room Size','Label','meters','Mapping',{'int',1,50}),...
        % 
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
            rIndexZ5 = wIndex - plugin.NSamples(5);
            rIndexZ6 = wIndex - plugin.NSamples(6);
            rIndexZ7 = wIndex - plugin.NSamples(7);
            rIndexZ8 = wIndex - plugin.NSamples(8);
            rIndexZ9 = wIndex - plugin.NSamples(9);
            rIndexZ10 = wIndex - plugin.NSamples(10);
            rIndexZ11 = wIndex - plugin.NSamples(11);
            rIndexZ12 = wIndex - plugin.NSamples(12);
            rIndexZ13 = wIndex - plugin.NSamples(13);
            rIndexZ14 = wIndex - plugin.NSamples(14);
            rIndexZ15 = wIndex - plugin.NSamples(15);
            rIndexZ16 = wIndex - plugin.NSamples(16);
            
            rIndexTDL_1 = wIndexTap - plugin.TSamples(1);
            rIndexTDL_2 = wIndexTap - plugin.TSamples(2);
            rIndexTDL_3 = wIndexTap - plugin.TSamples(3);
            rIndexTDL_4 = wIndexTap - plugin.TSamples(4);
            rIndexTDL_5 = wIndexTap - plugin.TSamples(5);
            rIndexTDL_6 = wIndexTap - plugin.TSamples(6);
            
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
            if rIndexZ5 <= 0 
                rIndexZ5 = rIndexZ5 + 192001;
            end
            if rIndexZ6 <= 0 
                rIndexZ6 = rIndexZ6 + 192001;
            end
            if rIndexZ7 <= 0 
                rIndexZ7 = rIndexZ7 + 192001;
            end
            if rIndexZ8 <= 0 
                rIndexZ8 = rIndexZ8 + 192001;
            end
            if rIndexZ9 <= 0 
                rIndexZ9 = rIndexZ9 + 192001;
            end
            if rIndexZ10 <= 0 
                rIndexZ10 = rIndexZ10 + 192001;
            end
            if rIndexZ11 <= 0 
                rIndexZ11 = rIndexZ11 + 192001;
            end
            if rIndexZ12 <= 0 
                rIndexZ12 = rIndexZ12 + 192001;
            end
            if rIndexZ13 <= 0 
                rIndexZ13 = rIndexZ13 + 192001;
            end
            if rIndexZ14 <= 0 
                rIndexZ14 = rIndexZ14 + 192001;
            end
            if rIndexZ15 <= 0 
                rIndexZ15 = rIndexZ15 + 192001;
            end
            if rIndexZ16 <= 0 
                rIndexZ16 = rIndexZ16 + 192001;
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
            if rIndexTDL_5 <= 0
                rIndexTDL_5 = rIndexTDL_5 + 192001;
            end
            if rIndexTDL_6 <= 0
                rIndexTDL_6 = rIndexTDL_6 + 192001;
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
                          KK(4)*plugin.tdl(rIndexTDL_4,:) + ...
                          KK(5)*plugin.tdl(rIndexTDL_5,:) + ...
                          KK(6)*plugin.tdl(rIndexTDL_6,:);
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
                        plugin.z3(rIndexZ3) plugin.z4(rIndexZ4)...
                        plugin.z5(rIndexZ5) plugin.z6(rIndexZ6)...
                        plugin.z7(rIndexZ7) plugin.z8(rIndexZ8)...
                        plugin.z9(rIndexZ9) plugin.z10(rIndexZ10)...
                        plugin.z11(rIndexZ11) plugin.z12(rIndexZ12)...
                        plugin.z13(rIndexZ13) plugin.z14(rIndexZ14)...
                        plugin.z15(rIndexZ15) plugin.z16(rIndexZ16)];    
                % ---------------------------------------------------------
                % LowPass Filters after each delay line
                plugin.lpfPrev(1) = plugin.B0*plugin.z1(rIndexZ1) + B1*plugin.lpfPrev(1);
                plugin.lpfPrev(2) = plugin.B0*plugin.z2(rIndexZ2) + B1*plugin.lpfPrev(2);
                plugin.lpfPrev(3) = plugin.B0*plugin.z3(rIndexZ3) + B1*plugin.lpfPrev(3);
                plugin.lpfPrev(4) = plugin.B0*plugin.z4(rIndexZ4) + B1*plugin.lpfPrev(4);
                plugin.lpfPrev(5) = plugin.B0*plugin.z5(rIndexZ5) + B1*plugin.lpfPrev(5);
                plugin.lpfPrev(6) = plugin.B0*plugin.z6(rIndexZ6) + B1*plugin.lpfPrev(6);
                plugin.lpfPrev(7) = plugin.B0*plugin.z7(rIndexZ7) + B1*plugin.lpfPrev(7);
                plugin.lpfPrev(8) = plugin.B0*plugin.z8(rIndexZ8) + B1*plugin.lpfPrev(8);
                plugin.lpfPrev(9) = plugin.B0*plugin.z9(rIndexZ9) + B1*plugin.lpfPrev(9);
                plugin.lpfPrev(10) = plugin.B0*plugin.z10(rIndexZ10) + B1*plugin.lpfPrev(10);
                plugin.lpfPrev(11) = plugin.B0*plugin.z11(rIndexZ11) + B1*plugin.lpfPrev(11);
                plugin.lpfPrev(12) = plugin.B0*plugin.z12(rIndexZ12) + B1*plugin.lpfPrev(12);
                plugin.lpfPrev(13) = plugin.B0*plugin.z13(rIndexZ13) + B1*plugin.lpfPrev(13);
                plugin.lpfPrev(14) = plugin.B0*plugin.z14(rIndexZ14) + B1*plugin.lpfPrev(14);
                plugin.lpfPrev(15) = plugin.B0*plugin.z15(rIndexZ15) + B1*plugin.lpfPrev(15);
                plugin.lpfPrev(16) = plugin.B0*plugin.z16(rIndexZ16) + B1*plugin.lpfPrev(16);
                % ---------------------------------------------------------
                % equation -- output
                out(i,1) = yl + CC(1)*plugin.lpfPrev(1) + CC(2)*plugin.lpfPrev(2) ...
                              + CC(3)*plugin.lpfPrev(3) + CC(4)*plugin.lpfPrev(4) ...
                              + CC(5)*plugin.lpfPrev(5) + CC(6)*plugin.lpfPrev(6) ...
                              + CC(7)*plugin.lpfPrev(7) + CC(8)*plugin.lpfPrev(8) ...
                              + CC(9)*plugin.lpfPrev(9) + CC(10)*plugin.lpfPrev(10) ...
                              + CC(11)*plugin.lpfPrev(11) + CC(12)*plugin.lpfPrev(12) ...
                              + CC(13)*plugin.lpfPrev(13) + CC(14)*plugin.lpfPrev(14) ...
                              + CC(15)*plugin.lpfPrev(15) + CC(16)*plugin.lpfPrev(16);
                        
                out(i,2) = yr + CC(1)*plugin.lpfPrev(1) + CC(2)*plugin.lpfPrev(2) ...
                              + CC(3)*plugin.lpfPrev(3) + CC(4)*plugin.lpfPrev(4) ...
                              + CC(5)*plugin.lpfPrev(5) + CC(6)*plugin.lpfPrev(6) ...
                              + CC(7)*plugin.lpfPrev(7) + CC(8)*plugin.lpfPrev(8) ...
                              + CC(9)*plugin.lpfPrev(9) + CC(10)*plugin.lpfPrev(10) ...
                              + CC(11)*plugin.lpfPrev(11) + CC(12)*plugin.lpfPrev(12) ...
                              + CC(13)*plugin.lpfPrev(13) + CC(14)*plugin.lpfPrev(14) ...
                              + CC(15)*plugin.lpfPrev(15) + CC(16)*plugin.lpfPrev(16);
                % ---------------------------------------------------------
                % -- write to circular buffer
                % plugin.tdl(wIndexTap,:) = in(i,:) * plugin.Gain;
                plugin.tdl(wIndexTap,:) = in(i,:) * (2/plugin.distance);
                plugin.bufferL(wIndexL) = in(i,1) * (2/plugin.distance) + tdl_out(1,1);
                plugin.bufferR(wIndexR) = in(i,2) * (2/plugin.distance) + tdl_out(1,2);
                
                plugin.z1(wIndex,:) = plugin.tdl(rIndexTDL_1,:)*BB(1) + temp*A(1,:)';
                plugin.z2(wIndex,:) = plugin.tdl(rIndexTDL_2,:)*BB(2) + temp*A(2,:)';
                plugin.z3(wIndex,:) = plugin.tdl(rIndexTDL_3,:)*BB(3) + temp*A(3,:)';
                plugin.z4(wIndex,:) = plugin.tdl(rIndexTDL_4,:)*BB(4) + temp*A(4,:)';
                plugin.z5(wIndex,:) = plugin.tdl(rIndexTDL_5,:)*BB(5) + temp*A(5,:)';
                plugin.z6(wIndex,:) = plugin.tdl(rIndexTDL_6,:)*BB(6) + temp*A(6,:)';
                plugin.z7(wIndex,:) = plugin.tdl(rIndexTDL_1,:)*BB(7) + temp*A(7,:)';
                plugin.z8(wIndex,:) = plugin.tdl(rIndexTDL_2,:)*BB(8) + temp*A(8,:)';
                plugin.z9(wIndex,:) = plugin.tdl(rIndexTDL_3,:)*BB(9) + temp*A(9,:)';
                plugin.z10(wIndex,:) = plugin.tdl(rIndexTDL_4,:)*BB(10) + temp*A(10,:)';
                plugin.z11(wIndex,:) = plugin.tdl(rIndexTDL_5,:)*BB(11) + temp*A(11,:)';
                plugin.z12(wIndex,:) = plugin.tdl(rIndexTDL_6,:)*BB(12) + temp*A(12,:)';
                plugin.z13(wIndex,:) = plugin.tdl(rIndexTDL_1,:)*BB(13) + temp*A(13,:)';
                plugin.z14(wIndex,:) = plugin.tdl(rIndexTDL_2,:)*BB(14) + temp*A(14,:)';
                plugin.z15(wIndex,:) = plugin.tdl(rIndexTDL_3,:)*BB(15) + temp*A(15,:)';
                plugin.z16(wIndex,:) = plugin.tdl(rIndexTDL_4,:)*BB(16) + temp*A(16,:)';
                
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
                rIndexZ5 = rIndexZ5 + 1;
                rIndexZ6 = rIndexZ6 + 1;
                rIndexZ7 = rIndexZ7 + 1;
                rIndexZ8 = rIndexZ8 + 1;
                rIndexZ9 = rIndexZ9 + 1;
                rIndexZ10 = rIndexZ10 + 1;
                rIndexZ11 = rIndexZ11 + 1;
                rIndexZ12 = rIndexZ12 + 1;
                rIndexZ13 = rIndexZ13 + 1;
                rIndexZ14 = rIndexZ14 + 1;
                rIndexZ15 = rIndexZ15 + 1;
                rIndexZ16 = rIndexZ16 + 1;
                
                rIndexTDL_1 = rIndexTDL_1 + 1;
                rIndexTDL_2 = rIndexTDL_2 + 1;
                rIndexTDL_3 = rIndexTDL_3 + 1;
                rIndexTDL_4 = rIndexTDL_4 + 1;
                rIndexTDL_5 = rIndexTDL_5 + 1;
                rIndexTDL_6 = rIndexTDL_6 + 1;
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
                if rIndexZ5 > 192001
                    rIndexZ5 = 1;
                end
                if rIndexZ6 > 192001
                    rIndexZ6 = 1;
                end
                if rIndexZ7 > 192001
                    rIndexZ7 = 1;
                end
                if rIndexZ8 > 192001
                    rIndexZ8 = 1;
                end
                if rIndexZ9 > 192001
                    rIndexZ9 = 1;
                end
                if rIndexZ10 > 192001
                    rIndexZ10 = 1;
                end
                if rIndexZ11 > 192001
                    rIndexZ11 = 1;
                end
                if rIndexZ12 > 192001
                    rIndexZ12 = 1;
                end
                if rIndexZ13 > 192001
                    rIndexZ13 = 1;
                end
                if rIndexZ14 > 192001
                    rIndexZ14 = 1;
                end
                if rIndexZ15 > 192001
                    rIndexZ15 = 1;
                end
                if rIndexZ16 > 192001
                    rIndexZ16 = 1;
                end
                % ---------------------------------------------------------
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
                if rIndexTDL_5 > 192001
                    rIndexTDL_5 = 1;
                end
                if rIndexTDL_6 > 192001
                    rIndexTDL_6 = 1;
                end
                % ---------------------------------------------------------
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
%         function set.pathmin(plugin, val)
%             fs = getSampleRate(plugin);
%             Np = 16;
%             i = [1:Np];
%             % Approximate desired delay-line lengths using powers of distinct primes:
%             % c = 343; % soundspeed in m/s at 20 degrees C for dry air
%             plugin.pathmin = val;
%             plugin.dmin = fs*val/plugin.c;
%             dl = plugin.dmin * (plugin.dmax/plugin.dmin).^(i/(Np-1)); % desired delay in samples
%             ppwr = floor(0.5 + log(dl)./log(plugin.prime(i))); % best prime power
%             plugin.NSamples(i) = plugin.prime(i).^ppwr; % each delay a power of a distinct prime
%             ptr = '--------------------------------------------------------'
%             plugin.NSamples(i)
%             plugin.TSamples(i)
%             plugin.NSamplesL
%             plugin.NSamplesR
%         end
        %------------------------------------------------------------------
        function set.pathmax(plugin, val)
            Np = 16;
            i = [1:Np];
            plugin.pathmax = val;
            plugin.dmin = getSampleRate(plugin)*5/plugin.c;
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
