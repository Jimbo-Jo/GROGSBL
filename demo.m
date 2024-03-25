clc;
clear;
close all
font_size = 15;
line_width = 1.5;

%% Parameters
array_num = 16;                  % number of sensors
snapshot_num = 50;               % Snapshots
source_doa = [77.5,90.3,100.7];  % directions of arrival
c = 1520;                        % speed of sound
f = 500;                         % frequency
lambda = c/f;                    % wavelength
d = 0.5*lambda;                  % intersensor spacing
source_num = length(source_doa); % number of signal
snr = 0;                         % SNR
reso = 8;                        % grid resolution
reso_grid = (0:reso:180);        % initial grid

%% Signal generate
Asignal = exp(-1i*(0:array_num-1)'*2*pi*(d/lambda)*cosd(source_doa));
Xsource = exp(1i*2*pi*rand(source_num,snapshot_num));
Ysignal = Asignal*Xsource; 
Y=awgn(Ysignal,snr,'measured');

%% GROGSBL
etc = 3;
tic
[Pm_grogsbl,theta_grogsbl]=GROGSBL(Y,d,lambda,reso_grid,etc);
Pm_grogsbl = Pm_grogsbl/max(Pm_grogsbl);
Pm_grogsbl = 10*log10(Pm_grogsbl);
toc

%% Plot
figure
plot(theta_grogsbl,Pm_grogsbl,'LineWidth',line_width);
hold on;
plot(source_doa,max(Pm_grogsbl),'ro','LineWidth',line_width);
set(gca,'XLim',[0,180],'YLim',[-60,0]);
legend('GROGSBL','True DOAs','fontsize',font_size,'fontname','Times New Roman');
xlabel('\theta (\circ)','fontsize',font_size,'fontname','Times New Roman');
ylabel('Normalized Spectrum(dB)','fontsize',font_size,'fontname','Times New Roman');
hold off;
