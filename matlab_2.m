%% Analysis of 2D state of stress 
clc
clear all
close all
%% Define stress components
sigma_xx = 100; % units : Mpa
sigma_yy = -40; % units : Mpa
tau_xy = 30; % units  : Mpa
%% Define angle of rotation
theta = 30; % units : degrees
%% Calculate stress in coordinate axes 
sigma_xx_prime = ((sigma_xx + sigma_yy)/2) + (((sigma_xx - sigma_yy)/2).*cosd(2.*theta)) + tau_xy.*sind(2.*theta);
sigma_yy_prime = ((sigma_xx + sigma_yy)/2) - (((sigma_xx - sigma_yy)/2).*cosd(2.*theta)) - tau_xy.*sind(2.*theta);
tau_xy_prime = ((sigma_yy - sigma_xx)/2).*sind(2.*theta) + tau_xy.*cosd(2.*theta);
%% check
 sigma = [sigma_xx tau_xy; tau_xy sigma_yy]
 Q = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]
 sigma_prime = Q*sigma*transpose(Q)
 %% Calculate principal stresses
 Average = (sigma_xx + sigma_yy)/2
 radius = sqrt((((sigma_xx - sigma_yy)/2)^2) + (tau_xy)^2)
 sigma_max = Average + radius
 sigma_min = Average - radius
 tau_max = radius
 %% check 
 sigma_normal = eig(sigma)
 %% Calculation of theta_n and theta_s
 theta_n_1 = 0.5*atand((2*tau_xy)/(sigma_xx - sigma_yy))
 theta_n_2 = theta_n_1 + 90
 theta_s_1 = 0.5*atand(-((sigma_xx - sigma_yy))/(2*tau_xy))
 theta_s_2 = theta_s_1 + 90
 %% Analysis of 2D state of stress 
clc
clear all
close all
%% Define stress components
sigma_xx = 100; % units : Mpa
sigma_yy = -40; % units : Mpa
tau_xy = 30; % units  : Mpa
%% Define angle of rotation
theta = 30; % units : degrees
%% Calculate stress in coordinate axes 
sigma_xx_prime = ((sigma_xx + sigma_yy)/2) + (((sigma_xx - sigma_yy)/2).*cosd(2.*theta)) + tau_xy.*sind(2.*theta);
sigma_yy_prime = ((sigma_xx + sigma_yy)/2) - (((sigma_xx - sigma_yy)/2).*cosd(2.*theta)) - tau_xy.*sind(2.*theta);
tau_xy_prime = ((sigma_yy - sigma_xx)/2).*sind(2.*theta) + tau_xy.*cosd(2.*theta);
%% check
 %Sigma = [sigma_xx tau_xy; tau_xy sigma_xy]
 %Q = [cosd(theta) sind(theta); -sind(theta) cosd(theta)]
 %sigma_prime = Q*sigma*transpose(Q)
  %% Calculate principal stresses
 Average = (sigma_xx + sigma_yy)/2
 radius = sqrt((((sigma_xx - sigma_yy)/2)^2) + (tau_xy)^2)
 sigma_max = Average + radius
 sigma_min = Average - radius
 tau_max = radius
  %% check 
 % sigma_normal = eig(Sigma)
  %% Calculation of theta_n and theta_s
 theta_n_1 = 0.5*atand((2*tau_xy)/(sigma_xx - sigma_yy))
 theta_n_2 = theta_n_1 + 90
 theta_s_1 = 0.5*atand(-((sigma_xx - sigma_yy))/(2*tau_xy))
 theta_s_2 = theta_s_1 + 90


 %% Plot
 figure(1)
 plot(theta,sigma_xx_prime,'r',theta,sigma_yy_prime,'b',theta,tau_xy_prime,'g')
 xlabel('Angle(degrees)')
 ylabel('stress (Mpa)')
 grid on 
 grid minor
% saveas(figure(1),'Plot_stress_vs_theta.jpg')
 legend('sigmaprime_{xx}','sigmaprime_{yy}','tauprime_{xy}')
 saveas(figure(1),'plot_Stress_vs_Theta.jpeg')
%% 
[sigma_xx_prime_max,M] = max(sigma_xx_prime)
theta_xx_max = theta(M+1)

