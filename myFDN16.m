%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Fri Mar 19 14:30:18 CET 2017
% Last Modified by  : Matteo Girardi (girardi.matthew@gmail.com)
% Last Modified on  : Wed Apr 26 23:57:09 CEST 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ~~~~~~~~~~~~~~~ -*- Feedback Delay Network -*- ~~~~~~~~~~~~~~~~~~~~~~ %%
% Real-time implementation of FDN
% you need the Audio System Toolbox!! 
% 
% 16 Delay lines + LOWPASS Filter
% based on: 
% - Physical Audio Signal Processing
%   for Virtual Musical Instruments and Audio Effects
%   Julius O. Smith III
% p. 65-67, p.85-127
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef myFDN16 < audioPlugin
    properties 
        % LPF Coeff
        B0 = 0.1;
        % Coeff before and after del_buffer
        B = rand(1,16);
        C = rand(1,16);
        bM = 0.5;
        cM = 0.5;
        % Feedback Matrix
        Dampening = 0.1;
        % Input Gain
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
        % tapped delay line coeff k, random
        k = rand(1,16);
        kERefl = 0.5;
    end
    properties (Access = private)
        % Tapped Delay line
        tapDel = zeros(44100,2);
        % Delay Lines
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
        % Index
        BufferIndex = 1;
        BufferTapIndex = 1;
        % Tapped Delay Lines 
        TSamples = [432 464 476 570 635 683 747 802 867 922 995 1048 1148 1170 1181 1192];
        % Delay times
        NSamples = [512 729 625 343 1331 2197 289 361 529 841 961 1369 1681 1849 2209 2809];
        % Last output
        lpfPrev = zeros(1,16);     
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Gain','DisplayName','Dry','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('cM','DisplayName','Wet','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('kERefl','DisplayName','Pre-reverb','Label','%','Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening','DisplayName','Dampening','Label','%','Mapping',{'lin',0,0.5}),...
            audioPluginParameter('B0','DisplayName','Lowpass','Mapping',{'lin',0,1}),...
            audioPluginParameter('pathmax','DisplayName','Room Size','Label','meters','Mapping',{'int',5,50}));
            % audioPluginParameter('pathmin','DisplayName','Min Room Size','Label','meters','Mapping',{'int',1,50}),...
    end
    methods
        %------------------------------------------------------------------
        function out = process(plugin, in)
            % 
            out = zeros(size(in));
            % 
            writeIndex = plugin.BufferIndex;
            writeTap = plugin.BufferTapIndex;
            % wrapping
            td1_readIndex = writeTap - plugin.TSamples(1);
            td2_readIndex = writeTap - plugin.TSamples(2);
            td3_readIndex = writeTap - plugin.TSamples(3);
            td4_readIndex = writeTap - plugin.TSamples(4);
            td5_readIndex = writeTap - plugin.TSamples(5);
            td6_readIndex = writeTap - plugin.TSamples(6);
            td7_readIndex = writeTap - plugin.TSamples(7);
            td8_readIndex = writeTap - plugin.TSamples(8);
            td9_readIndex = writeTap - plugin.TSamples(9);
            td10_readIndex = writeTap - plugin.TSamples(10);
            td11_readIndex = writeTap - plugin.TSamples(11);
            td12_readIndex = writeTap - plugin.TSamples(12);
            td13_readIndex = writeTap - plugin.TSamples(13);
            td14_readIndex = writeTap - plugin.TSamples(14);
            td15_readIndex = writeTap - plugin.TSamples(15);
            td16_readIndex = writeTap - plugin.TSamples(16);
            % wrapping 
            Z1_readIndex = writeIndex - plugin.NSamples(1);
            Z2_readIndex = writeIndex - plugin.NSamples(2);
            Z3_readIndex = writeIndex - plugin.NSamples(3);
            Z4_readIndex = writeIndex - plugin.NSamples(4);
            Z5_readIndex = writeIndex - plugin.NSamples(5);
            Z6_readIndex = writeIndex - plugin.NSamples(6);
            Z7_readIndex = writeIndex - plugin.NSamples(7);
            Z8_readIndex = writeIndex - plugin.NSamples(8);
            Z9_readIndex = writeIndex - plugin.NSamples(9);
            Z10_readIndex = writeIndex - plugin.NSamples(10);
            Z11_readIndex = writeIndex - plugin.NSamples(11);
            Z12_readIndex = writeIndex - plugin.NSamples(12);
            Z13_readIndex = writeIndex - plugin.NSamples(13);
            Z14_readIndex = writeIndex - plugin.NSamples(14);
            Z15_readIndex = writeIndex - plugin.NSamples(15);
            Z16_readIndex = writeIndex - plugin.NSamples(16);
            % ---------------------------------------------------------
            %
            if Z1_readIndex <= 0
                Z1_readIndex = Z1_readIndex + 192001;
            end
            if Z2_readIndex <= 0
                Z2_readIndex = Z2_readIndex + 192001;
            end
            if Z3_readIndex <= 0
                Z3_readIndex = Z3_readIndex + 192001;
            end
            if Z4_readIndex <= 0
                Z4_readIndex = Z4_readIndex + 192001;
            end
            if Z5_readIndex <= 0
                Z5_readIndex = Z5_readIndex + 192001;
            end
            if Z6_readIndex <= 0
                Z6_readIndex = Z6_readIndex + 192001;
            end
            if Z7_readIndex <= 0
                Z7_readIndex = Z7_readIndex + 192001;
            end
            if Z8_readIndex <= 0
                Z8_readIndex = Z8_readIndex + 192001;
            end
            if Z9_readIndex <= 0
                Z9_readIndex = Z9_readIndex + 192001;
            end
            if Z10_readIndex <= 0
                Z10_readIndex = Z10_readIndex + 192001;
            end
            if Z11_readIndex <= 0
                Z11_readIndex = Z11_readIndex + 192001;
            end
            if Z12_readIndex <= 0
                Z12_readIndex = Z12_readIndex + 192001;
            end
            if Z13_readIndex <= 0
                Z13_readIndex = Z13_readIndex + 192001;
            end
            if Z14_readIndex <= 0
                Z14_readIndex = Z14_readIndex + 192001;
            end
            if Z15_readIndex <= 0
                Z15_readIndex = Z15_readIndex + 192001;
            end
            if Z16_readIndex <= 0
                Z16_readIndex = Z16_readIndex + 192001;
            end
            % ---------------------------------------------------------
            % 
            if td1_readIndex <= 0
                td1_readIndex = td1_readIndex + 44100;
            end
            if td2_readIndex <= 0
                td2_readIndex = td2_readIndex + 44100;
            end
            if td3_readIndex <= 0
                td3_readIndex = td3_readIndex + 44100;
            end
            if td4_readIndex <= 0
                td4_readIndex = td4_readIndex + 44100;
            end
            if td5_readIndex <= 0
                td5_readIndex = td5_readIndex + 44100;
            end
            if td6_readIndex <= 0
                td6_readIndex = td6_readIndex + 44100;
            end
            if td7_readIndex <= 0
                td7_readIndex = td7_readIndex + 44100;
            end
            if td8_readIndex <= 0
                td8_readIndex = td8_readIndex + 44100;
            end
            if td9_readIndex <= 0
                td9_readIndex = td9_readIndex + 44100;
            end
            if td1_readIndex <= 0
                td1_readIndex = td1_readIndex + 44100;
            end
            if td10_readIndex <= 0
                td10_readIndex = td10_readIndex + 44100;
            end
            if td11_readIndex <= 0
                td11_readIndex = td11_readIndex + 44100;
            end
            if td12_readIndex <= 0
                td12_readIndex = td12_readIndex + 44100;
            end
            if td13_readIndex <= 0
                td13_readIndex = td13_readIndex + 44100;
            end
            if td14_readIndex <= 0
                td14_readIndex = td14_readIndex + 44100;
            end
            if td15_readIndex <= 0
                td15_readIndex = td15_readIndex + 44100;
            end
            if td16_readIndex <= 0
                td16_readIndex = td16_readIndex + 44100;
            end
            
            for i = 1:size(in,1)
                % b and c coff - buffers
                bN = plugin.B*plugin.bM;
                cN = plugin.C*plugin.cM;
                % k coeff, tapped del lines
                kN = plugin.k*plugin.kERefl;
                % LPF coeff
                B1 = 1 - plugin.B0;
                % feedback matrix
                A = plugin.Dampening*(1/2)*hadamard(16);
                % -------------------------------------
                % tapped delay line
                tap_out = kN(1)*plugin.tapDel(td1_readIndex) + ...
                          kN(2)*plugin.tapDel(td2_readIndex) + ...
                          kN(3)*plugin.tapDel(td3_readIndex) + ...
                          kN(4)*plugin.tapDel(td4_readIndex) + ...
                          kN(5)*plugin.tapDel(td5_readIndex) + ...
                          kN(6)*plugin.tapDel(td6_readIndex) + ...
                          kN(7)*plugin.tapDel(td7_readIndex) + ...
                          kN(8)*plugin.tapDel(td8_readIndex) + ...
                          kN(9)*plugin.tapDel(td9_readIndex) + ...
                          kN(10)*plugin.tapDel(td10_readIndex) + ...
                          kN(11)*plugin.tapDel(td11_readIndex) + ...
                          kN(12)*plugin.tapDel(td12_readIndex) + ...
                          kN(13)*plugin.tapDel(td13_readIndex) + ...
                          kN(14)*plugin.tapDel(td14_readIndex) + ...
                          kN(15)*plugin.tapDel(td15_readIndex) + ...
                          kN(16)*plugin.tapDel(td16_readIndex);
                % -------------------------------------
                % tmp delay samples
                temp = [plugin.z1(Z1_readIndex)...
                        plugin.z2(Z2_readIndex)...
                        plugin.z3(Z3_readIndex)...
                        plugin.z4(Z4_readIndex)...
                        plugin.z5(Z4_readIndex)...
                        plugin.z6(Z6_readIndex)...
                        plugin.z7(Z7_readIndex)...
                        plugin.z8(Z8_readIndex)...
                        plugin.z9(Z9_readIndex)...
                        plugin.z10(Z10_readIndex)...
                        plugin.z11(Z11_readIndex)...
                        plugin.z12(Z12_readIndex)...
                        plugin.z13(Z13_readIndex)...
                        plugin.z14(Z14_readIndex)...
                        plugin.z15(Z15_readIndex)...
                        plugin.z16(Z16_readIndex)];
                % -------------------------------------
                % LowPass Filters after each delay line
                plugin.lpfPrev(1) = plugin.B0*plugin.z1(Z1_readIndex) + B1*plugin.lpfPrev(1);
                plugin.lpfPrev(2) = plugin.B0*plugin.z2(Z2_readIndex) + B1*plugin.lpfPrev(2);
                plugin.lpfPrev(3) = plugin.B0*plugin.z3(Z3_readIndex) + B1*plugin.lpfPrev(3);
                plugin.lpfPrev(4) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(4);
                plugin.lpfPrev(5) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(5);
                plugin.lpfPrev(6) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(6);
                plugin.lpfPrev(7) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(7);
                plugin.lpfPrev(8) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(8);
                plugin.lpfPrev(9) = plugin.B0*plugin.z1(Z1_readIndex) + B1*plugin.lpfPrev(9);
                plugin.lpfPrev(10) = plugin.B0*plugin.z2(Z2_readIndex) + B1*plugin.lpfPrev(10);
                plugin.lpfPrev(11) = plugin.B0*plugin.z3(Z3_readIndex) + B1*plugin.lpfPrev(11);
                plugin.lpfPrev(12) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(12);
                plugin.lpfPrev(13) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(13);
                plugin.lpfPrev(14) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(14);
                plugin.lpfPrev(15) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(15);
                plugin.lpfPrev(16) = plugin.B0*plugin.z4(Z4_readIndex) + B1*plugin.lpfPrev(16);
                % -------------------------------------
                % equation -- output
                out(i,:) = (in(i,:)*plugin.Gain) + tap_out + ...
                            cN(1)*plugin.lpfPrev(1) + ...
                            cN(2)*plugin.lpfPrev(2) + ...
                            cN(3)*plugin.lpfPrev(3) + ...
                            cN(4)*plugin.lpfPrev(4) + ...
                            cN(5)*plugin.lpfPrev(5) + ...
                            cN(6)*plugin.lpfPrev(6) + ...
                            cN(7)*plugin.lpfPrev(7) + ...
                            cN(8)*plugin.lpfPrev(8) + ...
                            cN(9)*plugin.lpfPrev(9) + ...
                            cN(10)*plugin.lpfPrev(10) + ...
                            cN(11)*plugin.lpfPrev(11) + ...
                            cN(12)*plugin.lpfPrev(12) + ...
                            cN(13)*plugin.lpfPrev(13) + ...
                            cN(14)*plugin.lpfPrev(14) + ...
                            cN(15)*plugin.lpfPrev(15) + ...
                            cN(16)*plugin.lpfPrev(16);
                % -------------------------------------
                % wirte to circular buffers
                plugin.tapDel(writeTap,:) = in(i,:) * plugin.Gain;

                plugin.z1(writeIndex,:) = in(i,:)*bN(1) + temp*A(1,:)';
                plugin.z2(writeIndex,:) = in(i,:)*bN(2) + temp*A(2,:)';
                plugin.z3(writeIndex,:) = in(i,:)*bN(3) + temp*A(3,:)';
                plugin.z4(writeIndex,:) = in(i,:)*bN(4) + temp*A(4,:)';
                plugin.z5(writeIndex,:) = in(i,:)*bN(5) + temp*A(5,:)';
                plugin.z6(writeIndex,:) = in(i,:)*bN(6) + temp*A(6,:)';
                plugin.z7(writeIndex,:) = in(i,:)*bN(7) + temp*A(7,:)';
                plugin.z8(writeIndex,:) = in(i,:)*bN(8) + temp*A(8,:)';
                plugin.z9(writeIndex,:) = in(i,:)*bN(9) + temp*A(9,:)';
                plugin.z10(writeIndex,:) = in(i,:)*bN(10) + temp*A(10,:)';
                plugin.z11(writeIndex,:) = in(i,:)*bN(11) + temp*A(11,:)';
                plugin.z12(writeIndex,:) = in(i,:)*bN(12) + temp*A(12,:)';
                plugin.z13(writeIndex,:) = in(i,:)*bN(13) + temp*A(13,:)';
                plugin.z14(writeIndex,:) = in(i,:)*bN(14) + temp*A(14,:)';
                plugin.z15(writeIndex,:) = in(i,:)*bN(15) + temp*A(15,:)';
                plugin.z16(writeIndex,:) = in(i,:)*bN(16) + temp*A(16,:)';

                % it works, however too much feedback
%                 plugin.z1(writeIndex,:) = in(i,:)*bN(1) + plugin.lpfPrev*A(1,:)';
%                 plugin.z2(writeIndex,:) = in(i,:)*bN(2) + plugin.lpfPrev*A(2,:)';
%                 plugin.z3(writeIndex,:) = in(i,:)*bN(3) + plugin.lpfPrev*A(3,:)';
%                 plugin.z4(writeIndex,:) = in(i,:)*bN(4) + plugin.lpfPrev*A(4,:)';
%                 
%                 plugin.z5(writeIndex,:) = in(i,:)*bN(5) + plugin.lpfPrev*A(5,:)';
%                 plugin.z6(writeIndex,:) = in(i,:)*bN(6) + plugin.lpfPrev*A(6,:)';
%                 plugin.z7(writeIndex,:) = in(i,:)*bN(7) + plugin.lpfPrev*A(7,:)';
%                 plugin.z8(writeIndex,:) = in(i,:)*bN(8) + plugin.lpfPrev*A(8,:)';
%                 
%                 plugin.z9(writeIndex,:) = in(i,:)*bN(9) + plugin.lpfPrev*A(9,:)';
%                 plugin.z10(writeIndex,:) = in(i,:)*bN(10) + plugin.lpfPrev*A(10,:)';
%                 plugin.z11(writeIndex,:) = in(i,:)*bN(11) + plugin.lpfPrev*A(11,:)';
%                 plugin.z12(writeIndex,:) = in(i,:)*bN(12) + plugin.lpfPrev*A(12,:)';
%                 
%                 plugin.z13(writeIndex,:) = in(i,:)*bN(13) + plugin.lpfPrev*A(13,:)';
%                 plugin.z14(writeIndex,:) = in(i,:)*bN(14) + plugin.lpfPrev*A(14,:)';
%                 plugin.z15(writeIndex,:) = in(i,:)*bN(15) + plugin.lpfPrev*A(15,:)';
%                 plugin.z16(writeIndex,:) = in(i,:)*bN(16) + plugin.lpfPrev*A(16,:)';
                
                writeIndex = writeIndex + 1;
                writeTap = writeTap + 1;
               
                if writeIndex > 192001
                    writeIndex = 1;
                end
                if writeTap > 44100
                    writeTap = 1;
                end
                % ---------------------------------------------------------
                % update buffers read 
                Z1_readIndex = Z1_readIndex + 1;
                Z2_readIndex = Z2_readIndex + 1;
                Z3_readIndex = Z3_readIndex + 1;
                Z4_readIndex = Z4_readIndex + 1;
                Z5_readIndex = Z5_readIndex + 1;
                Z6_readIndex = Z6_readIndex + 1;
                Z7_readIndex = Z7_readIndex + 1;
                Z8_readIndex = Z8_readIndex + 1;
                Z9_readIndex = Z9_readIndex + 1;
                Z10_readIndex = Z10_readIndex + 1;
                Z11_readIndex = Z11_readIndex + 1;
                Z12_readIndex = Z12_readIndex + 1;
                Z13_readIndex = Z13_readIndex + 1;
                Z14_readIndex = Z14_readIndex + 1;
                Z15_readIndex = Z15_readIndex + 1;
                Z16_readIndex = Z16_readIndex + 1;
                % ---------------------------------------------------------
                % update tapped del line read 
                td1_readIndex = td1_readIndex + 1;
                td2_readIndex = td2_readIndex + 1;
                td3_readIndex = td3_readIndex + 1;
                td4_readIndex = td4_readIndex + 1;
                td5_readIndex = td5_readIndex + 1;
                td6_readIndex = td6_readIndex + 1;
                td7_readIndex = td7_readIndex + 1;
                td8_readIndex = td8_readIndex + 1;
                td9_readIndex = td9_readIndex + 1;
                td10_readIndex = td10_readIndex + 1;
                td11_readIndex = td11_readIndex + 1;
                td12_readIndex = td12_readIndex + 1;
                td13_readIndex = td13_readIndex + 1;
                td14_readIndex = td14_readIndex + 1;
                td15_readIndex = td15_readIndex + 1;
                td16_readIndex = td16_readIndex + 1;
                % ---------------------------------------------------------
                % delay lines
                if Z1_readIndex > 192001
                    Z1_readIndex = 1;
                end
                if Z2_readIndex > 192001
                    Z2_readIndex = 1;
                end
                if Z3_readIndex > 192001
                    Z3_readIndex = 1;
                end
                if Z4_readIndex > 192001
                    Z4_readIndex = 1;
                end
                if Z5_readIndex > 192001
                    Z5_readIndex = 1;
                end
                if Z6_readIndex > 192001
                    Z6_readIndex = 1;
                end
                if Z7_readIndex > 192001
                    Z7_readIndex = 1;
                end
                if Z8_readIndex > 192001
                    Z8_readIndex = 1;
                end
                if Z9_readIndex > 192001
                    Z9_readIndex = 1;
                end
                if Z10_readIndex > 192001
                    Z10_readIndex = 1;
                end
                if Z11_readIndex > 192001
                    Z11_readIndex = 1;
                end
                if Z12_readIndex > 192001
                    Z12_readIndex = 1;
                end
                if Z13_readIndex > 192001
                    Z13_readIndex = 1;
                end
                if Z14_readIndex > 192001
                    Z14_readIndex = 1;
                end
                if Z15_readIndex > 192001
                    Z15_readIndex = 1;
                end
                if Z16_readIndex > 192001
                    Z16_readIndex = 1;
                end
                % ---------------------------------------------------------
                % tapped Delay lines
                if td1_readIndex > 44100
                    td1_readIndex = 1;
                end
                if td2_readIndex > 44100
                    td2_readIndex = 1;
                end
                if td3_readIndex > 44100
                    td3_readIndex = 1;
                end
                if td4_readIndex > 44100
                    td4_readIndex = 1;
                end
                if td5_readIndex > 44100
                    td5_readIndex = 1;
                end
                if td6_readIndex > 44100
                    td6_readIndex = 1;
                end
                if td7_readIndex > 44100
                    td7_readIndex = 1;
                end
                if td8_readIndex > 44100
                    td8_readIndex = 1;
                end
                if td9_readIndex > 44100
                    td9_readIndex = 1;
                end
                if td10_readIndex > 44100
                    td10_readIndex = 1;
                end
                if td11_readIndex > 44100
                    td11_readIndex = 1;
                end
                if td12_readIndex > 44100
                    td12_readIndex = 1;
                end
                if td13_readIndex > 44100
                    td13_readIndex = 1;
                end
                if td14_readIndex > 44100
                    td14_readIndex = 1;
                end
                if td15_readIndex > 44100
                    td15_readIndex = 1;
                end
                if td16_readIndex > 44100
                    td16_readIndex = 1;
                end
            end
            plugin.BufferIndex = writeIndex;
            plugin.BufferTapIndex = writeTap;
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
%             ppwr = floor(log(dl)./log(plugin.prime(i))); % best prime power
%             plugin.NSamples = plugin.prime(i).^ppwr; % each delay a power of a distinct prime
%             plugin.NSamples
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
            plugin.NSamples = plugin.prime(i).^ppwr;
        end
        function reset(plugin)
            plugin.z1 = zeros(192001,2);
            plugin.z2 = zeros(192001,2); 
            plugin.z3 = zeros(192001,2); 
            plugin.z4 = zeros(192001,2); 
            plugin.z5 = zeros(192001,2); 
            plugin.z6 = zeros(192001,2); 
            plugin.z7 = zeros(192001,2); 
            plugin.z8 = zeros(192001,2); 
            plugin.z9 = zeros(192001,2); 
            plugin.z10 = zeros(192001,2); 
            plugin.z11 = zeros(192001,2); 
            plugin.z12 = zeros(192001,2); 
            plugin.z13 = zeros(192001,2); 
            plugin.z14 = zeros(192001,2); 
            plugin.z15 = zeros(192001,2); 
            plugin.z16 = zeros(192001,2); 
            plugin.tapDel = zeros(44100,2);
            plugin.NSamples = [512 729 625 343 1331 2197 289 361 529 841 961 1369 1681 1849 2209 2809];
            plugin.TSamples = [432 464 476 570 635 683 747 802 867 922 995 1048 1148 1170 1181 1192];
        end
    end
end
