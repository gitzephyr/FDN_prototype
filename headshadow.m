%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author            : Matteo Girardi
% Created on        : Sat Apr 15 12:40:26 CEST 2017
% Last Modified by  : Matteo Girardi (girardi dot matthew at gmail.com)
% Last Modified on  : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Head shadow model by Brown and Duda
% An Efficient Hrtf Model For 3-D Sound (1997)
% by C. Phillip Brown , Richard O. Duda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read soundfile and plot it
% SOUND FILE MUS BE MONO, BETTER IF RECORDED IN AN ANECHOIC CHAMBER
clear all; close all;
[buffer, fs] = audioread('snd/singing.wav');
dt = 1/fs;
tt = 0:dt:(length(buffer)*dt)-dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot it
plot(tt,buffer); xlabel('Seconds'); ylabel('Amplitude');
title('Opera voice');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% play it!
sound(buffer, fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Standard head radius 
std_circ = 58.5; % male, cm. From wikipedia
% std_circ = 53; % female, cm
a = std_circ/(2*pi); % cm
a = a/100; % meters
% 0.0875m (8.75 cm) is the standard head radius. > no reference
% a = 0.0875    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITD -- Lord Rayleigh
itd = [];
% a source in front of the head
% steps = [-90:5:90];
step = 0.001;
front_ang = [0:step:pi/2];
c = 343;
for theta = front_ang
    itd = [itd, (a/c)*(theta+sin(theta))];
end
% The model used had a distance between 
% the 2 ears of approximately 22?23 cm. 
% Initial measurements found that there was a maximum time delay
% of approximately 660 ?s (~0.7 ms) when the sound source was
% placed at directly 90° azimuth to one ear.
(a*2)/c;
0.23/c;
% a source at the back
back_ang = [pi/2:step:pi];
for theta = back_ang
    itd = [itd, (a/c)*(pi - theta+sin(theta))];
end
rad = [front_ang, back_ang];
plot(rad*(180/pi),itd); grid on;
title ('ITD - Interaural Time Difference');
xlabel('Theta (Degree°)'); ylabel('ITD (sec)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The Head model [Duda, Brown, 1997] 
deg_theta = [-180:10:180];
tl = [];
tr = [];
rad_theta = abs(deg_theta*(pi/180));
for i = rad_theta
    if i <= 1.5707
        tl = [tl, (a-a*sin(i))/c];
        tr = [tr, (a+a*(i))/c];
    else 
        tl = [tl, (a + a*(i*-1+pi))/c];
        tr = [tr, (a - a*sin(i)*-1)/c];
    end
end
plot(deg_theta,tl); hold on;
plot(deg_theta,tr); hold on;
title ('ITD - Interaural Time Difference');
xlabel('Theta (Degree°)'); ylabel('ITD (sec)');
legend('tl','tr'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING -- different alpha > start on the right
dT = [-90:20:90]; % theta 
% 0.0875m (8.75 cm) is the standard head radius,
beta = (2*343)/a; 
t = 1/fs;
tbeta = (t*beta);
for i = dT*(pi/180)
    alpha_l=1-cos(i+pi/2);
    alpha_r=1+cos(i+pi/2);

    % filter coefficients
    b0=2+tbeta; % b0 and b1 is unidndependend from theta and thus the same for both ears 
    b1=-2+tbeta;
    a0_l=2*alpha_l-tbeta; % a0 and a1 is dependend and thus different
    a1_l=-2*alpha_l+tbeta;
    a0_r=2*alpha_r-tbeta;
    a1_r=-2*alpha_r+tbeta;
    
    figure(1)
    freqz([b0,b1], [a0_l,a1_l]); hold on
    yl = filter([b0,b1], [a0_l,a1_l], buffer);
    title('Left ear');
    
    figure(2)
    freqz([b0,b1], [a0_r,a1_r]); hold on
    yr = filter([b0,b1], [a0_r,a1_r], buffer);
    title('Right ear');
    
    % soundsc([yr,yl],fs);
    % pause(5);

    figure(3);
    [h,w] = freqz([b0,b1], [a0_l,a1_l]);
    plot((abs(w)),20*log10(h)); hold on
    title('Left ear');

    figure(4);
    [h,w] = freqz([b0,b1], [a0_r,a1_r]);
    plot((abs(w)),20*log10(h)); hold on
    title('Right ear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TESTING -- alpha > start on the left
dT = [-90:20:90]; % theta 
% 0.0875m (8.75 cm) is the standard head radius,
beta = (2*343)/a; 
t = 1/fs;
tbeta = (t*beta);
for i = dT*(pi/180)
    alpha_l=1-sin(i);
    alpha_r=1+sin(i);

    % filter coefficients
    b0=2+tbeta; % b0 and b1 is unidndependend from theta and thus the same for both ears 
    b1=-2+tbeta;
    a0_l=2*alpha_l-tbeta; % a0 and a1 is dependend and thus different
    a1_l=-2*alpha_l+tbeta;
    a0_r=2*alpha_r-tbeta;
    a1_r=-2*alpha_r+tbeta;
    
    figure(1)
    freqz([b0,b1], [a0_l,a1_l]); hold on
    yl = filter([b0,b1], [a0_l,a1_l], buffer);
    title('Left ear');
    
    figure(2)
    freqz([b0,b1], [a0_r,a1_r]); hold on
    yr = filter([b0,b1], [a0_r,a1_r], buffer);
    title('Right ear');
    
    % soundsc([yr,yl],fs);
    % pause(5);

    figure(3);
    [h,w] = freqz([b0,b1], [a0_l,a1_l]);
    plot((abs(w)),20*log10(h)); hold on
    title('Left ear');

    figure(4);
    [h,w] = freqz([b0,b1], [a0_r,a1_r]);
    plot((abs(w)),20*log10(h)); hold on
    title('Right ear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% testing the filter > difference equation
% this is just in case of a stereofile
bufferl = buffer(:,1)'; 
l = length(bufferl);
% initialise
n = 1;
deg_theta = [-90:20:90];
dT = deg_theta;
for theta = dT
    yl=zeros(1,l);
    yr=zeros(1,l);
    alpha_l=1-sin(theta*(pi/180));
    alpha_r=1+sin(theta*(pi/180));

    % filter coefficients
    b0=2+tbeta; 
    b1=-2+tbeta;
    a0_l=2*alpha_l-tbeta;
    a1_l=-2*alpha_l+tbeta;
    a0_r=2*alpha_r-tbeta;
    a1_r=-2*alpha_r+tbeta;
    
    for i = 2:l
        yl(i)=((a0_l*bufferl(i))+(a1_l*bufferl(i-1)))-(b1*yl(i-1))*(1/b0);
        yr(i)=((a0_r*bufferl(i))+(a1_r*bufferl(i-1)))-(b1*yr(i-1))*(1/b0);
    end    
    % apply delay
    yl = [zeros(1,floor(tl(n)*fs)), yl];
    yr = [zeros(1,floor(tr(n)*fs)), yr];
    if length(yl) > length(yr)
        diff = length(yl) - length(yr);
        yr = [zeros(1,diff), yr];
    else
        diff = length(yr) - length(yl);
        yl = [zeros(1,diff), yl];
    end
    % plot it!
    figure(1)
    
    subplot(2,1,1);
    plot(yl);
    legend('Left ear');
    ylim([-2 2]);
    title(['Del Left: ', num2str(tl(n)*1000), ' ms. -- ', num2str(floor(tl(n)*fs)), 'samples']);
    ylabel('Amplitude'); xlabel('Time (sec)');
    
    subplot(2,1,2);
    plot(yr);
    title(['Del Rigth: ', num2str(tr(n)*1000), ' ms. -- ', num2str(floor(tr(n)*fs)), 'samples']);
    ylabel('Amplitude'); xlabel('Time (sec)');
    legend('Right ear');
    ylim([-2 2]);
    
    descr = {'Degree: ',num2str(theta), '°'};
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    ax2 = axes('Position',[1 1 1 1]);
    axes(ax1) % sets ax1 to current axes
    text(.025,0.6,descr)
    
    % play it!
    yout=[yl;yr];
    soundsc(yout,fs);
    pause(5);
    n = n + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%