classdef myFDN16 < audioPlugin
    properties 
        % LPF Coeff
        B0 = 0.97;
        % Delay Time
        Delay = 0.5;
        
        % Coeff before and after del_buffer
        B = 0.1;
        C = 0.1;
        % Feedback Matrix
        Dampening = 0.1;
        % Input Gain
        Gain = 0.1;
        % LPF
        yLast = 0;
    end
    properties (Access = private)
        % Delay Lines
        z1 = zeros(44100,2);
        z2 = zeros(44100,2);
        z3 = zeros(44100,2);
        z4 = zeros(44100,2);
        z5 = zeros(44100,2);
        z6 = zeros(44100,2);
        z7 = zeros(44100,2);
        z8 = zeros(44100,2);
        z9 = zeros(44100,2);
        z10 = zeros(44100,2);
        z11 = zeros(44100,2);
        z12 = zeros(44100,2);
        z13 = zeros(44100,2);
        z14 = zeros(44100,2);
        z15 = zeros(44100,2);
        z16 = zeros(44100,2);
        % 
        BufferIndex = 1;
        % Delay times
        % NSamples = [443 1949 4409 5417 6421 7537 8863 9049 10799 11177 12791 13679 14891 15287 16339 17657]';
        % NSamples = [32 243 625 343 1331 2197 4913 6859 12167 841 961 1369 1681 1849 2209 2809]';
        NSamples = [256 243 625 343 1331 2197 4913 6859 12167 841 961 1369 1681 1849 2209 2809]';
    end
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Gain',...
            'DisplayName','DRY',...
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('Dampening',...
            'DisplayName','Dampening',...
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('B',...
            'DisplayName','B Coeff',...
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('C',...
            'DisplayName','WET',...
            'Mapping',{'lin',0,1}),...
            audioPluginParameter('B0',...
            'DisplayName','LPF Coeff',...
            'Mapping',{'lin',0,1}));
%         PluginInterface = audioPluginInterface(...
%             audioPluginParameter('Delay',...
%             'DisplayName','Delay',...
%             'Label','seconds',...
%             'Mapping',{'lin',0,5}))

    end
    methods
        function out = process(plugin, in)
            % b and c coff - buffers
            bN = plugin.B*ones(1,16);
            cN = plugin.C*ones(1,16);
            
%             B1 = 1 - plugin.B0;
            
            out = zeros(size(in));
            writeIndex = plugin.BufferIndex;
            
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
            
            if Z1_readIndex <= 0
                Z1_readIndex = Z1_readIndex + 44100;
            end
            if Z2_readIndex <= 0
                Z2_readIndex = Z2_readIndex + 44100;
            end
            if Z3_readIndex <= 0
                Z3_readIndex = Z3_readIndex + 44100;
            end
            if Z4_readIndex <= 0
                Z4_readIndex = Z4_readIndex + 44100;
            end
            if Z5_readIndex <= 0
                Z5_readIndex = Z5_readIndex + 44100;
            end
            if Z6_readIndex <= 0
                Z6_readIndex = Z6_readIndex + 44100;
            end
            if Z7_readIndex <= 0
                Z7_readIndex = Z7_readIndex + 44100;
            end
            if Z8_readIndex <= 0
                Z8_readIndex = Z8_readIndex + 44100;
            end
            if Z9_readIndex <= 0
                Z9_readIndex = Z9_readIndex + 44100;
            end
            if Z10_readIndex <= 0
                Z10_readIndex = Z10_readIndex + 44100;
            end
            if Z11_readIndex <= 0
                Z11_readIndex = Z11_readIndex + 44100;
            end
            if Z12_readIndex <= 0
                Z12_readIndex = Z12_readIndex + 44100;
            end
            if Z13_readIndex <= 0
                Z13_readIndex = Z13_readIndex + 44100;
            end
            if Z14_readIndex <= 0
                Z14_readIndex = Z14_readIndex + 44100;
            end
            if Z15_readIndex <= 0
                Z15_readIndex = Z15_readIndex + 44100;
            end
            if Z16_readIndex <= 0
                Z16_readIndex = Z16_readIndex + 44100;
            end
            
            for i = 1:size(in,1)
                % feedback matrix
                A = plugin.Dampening*(1/2)*hadamard(16);
                
                temp = [plugin.z1(Z1_readIndex) plugin.z2(Z2_readIndex)...
                    plugin.z3(Z3_readIndex) plugin.z4(Z4_readIndex)...
                    plugin.z5(Z4_readIndex) plugin.z6(Z6_readIndex)...
                    plugin.z7(Z7_readIndex) plugin.z8(Z8_readIndex)...
                    plugin.z9(Z9_readIndex) plugin.z10(Z10_readIndex)...
                    plugin.z11(Z11_readIndex) plugin.z12(Z12_readIndex)...
                    plugin.z13(Z13_readIndex) plugin.z14(Z14_readIndex)...
                    plugin.z15(Z15_readIndex) plugin.z16(Z16_readIndex)];
                % equation
                y = (in(i,:) * plugin.Gain) + cN(1)*plugin.z1(Z1_readIndex) + ...
                    cN(2)*plugin.z2(Z2_readIndex) + cN(3)*plugin.z3(Z3_readIndex) + ...
                    cN(4)*plugin.z4(Z4_readIndex) + cN(5)*plugin.z5(Z5_readIndex) + ...
                    cN(6)*plugin.z6(Z6_readIndex) + cN(7)*plugin.z7(Z7_readIndex) + ...
                    cN(8)*plugin.z8(Z8_readIndex) + cN(9)*plugin.z9(Z9_readIndex) + ...
                    cN(10)*plugin.z10(Z10_readIndex) + cN(11)*plugin.z11(Z11_readIndex) + ...
                    cN(12)*plugin.z12(Z12_readIndex) + cN(13)*plugin.z13(Z13_readIndex) + ...
                    cN(14)*plugin.z14(Z14_readIndex) + cN(15)*plugin.z15(Z15_readIndex) + ...
                    cN(16)*plugin.z16(Z16_readIndex);
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
                
                writeIndex = writeIndex + 1;
                
                if writeIndex > 44100
                    writeIndex = 1;
                end
                
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
                
                if Z1_readIndex > 44100
                    Z1_readIndex = 1;
                end
                if Z2_readIndex > 44100
                    Z2_readIndex = 1;
                end
                if Z3_readIndex > 44100
                    Z3_readIndex = 1;
                end
                if Z4_readIndex > 44100
                    Z4_readIndex = 1;
                end
                
                if Z5_readIndex > 44100
                    Z5_readIndex = 1;
                end
                if Z6_readIndex > 44100
                    Z6_readIndex = 1;
                end
                if Z7_readIndex > 220500
                    Z7_readIndex = 1;
                end
                if Z8_readIndex > 44100
                    Z8_readIndex = 1;
                end
                
                if Z9_readIndex > 44100
                    Z9_readIndex = 1;
                end
                if Z10_readIndex > 44100
                    Z10_readIndex = 1;
                end
                if Z11_readIndex > 44100
                    Z11_readIndex = 1;
                end
                if Z12_readIndex > 44100
                    Z12_readIndex = 1;
                end
                if Z13_readIndex > 44100
                    Z13_readIndex = 1;
                end
                if Z14_readIndex > 44100
                    Z14_readIndex = 1;
                end
                if Z15_readIndex > 44100
                    Z15_readIndex = 1;
                end
                if Z16_readIndex > 44100
                    Z16_readIndex = 1;
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