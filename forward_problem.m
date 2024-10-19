clear all;
close all;
clc;

%% Parameters

% frequency of incident wave
f = 5*10^9;
c_0 = 3 * 10^8;
E_0 = 8.854 * 10^(-12);
M_0 = 4 * pi * 10^(-7);

% polarization
polarization = 'TM';

% first medium - vacuum
E1_r = 1;
M1_r = 1; 
n1 = sqrt(E1_r * M1_r);
E1 = E_0 * E1_r; 
k1 = 2*pi*f*n1/c_0;

% second medium - layer
E2_r = 1;
M2_r = 5; 
h = 1 * 10^-2;  % from cm to m
n2 = sqrt(E2_r * M2_r);
E2 = E_0 * E2_r;
k2 = 2*pi*f*n2/c_0;

% third medium - substrate
E3_r = 1.25;
M3_r = 1; 
n3 = sqrt(E3_r * M3_r);
E3 = E_0 * E3_r;
k3 = 2*pi*f*n3/c_0;

% define the ratio between first and third medium
ratio = sqrt((M3_r/E3_r)/(M1_r/E1_r));


%% Compute solution
counter = 0;

for theta = 0 : 90

    counter = counter+1;
    k1z = k1*cosd(theta);
    k2z = sqrt(k2^2 - (k1 * sind(theta))^2);
    k3z = sqrt(k3^2 - (k1 * sind(theta))^2);

    % solution in Matrix form
    A  = [1,          0,                       -1,                 -1;  ...
        0,      -exp(-1i*k3z*h),       exp(1i*k2z*h),         exp(- 1i*k2z*h); ...
        k1z/E1,       0,                      -k2z/E2,              k2z/E2; ...
        0,     (k3z/E3)*exp(-1i*k3z*h),  (k2z/E2)*exp(1i*k2z*h), (-k2z/E2)*exp(-1i*k2z*h)];

    C = [-1; 0; k1z/E1; 0];

    % compute solution for each incident angle -> assume Hi=1
    A_inv = inv(A);
    B = A_inv * C;   
    Bs(counter, :) = B;
  
    S1z_module(counter) =  (1/sqrt(M1_r/E1_r))*(1 - (abs(Bs(counter,1))).^2).* cosd(theta);
    S3z_module(counter) =  (1/sqrt(M3_r/E3_r))*ratio.^2*(abs(Bs(counter,2))).^2.* cos(asin(k1*sind(theta)/k3));

end


%% Plot transmission and reflection coefficients

theta= 0:90;

figure(1); 
grid on
hold on
reflection = plot(theta, abs(Bs(:,1)),'LineWidth', 2);
transmission = plot(theta, ratio*abs(Bs(:,2)),'LineWidth', 2);


legend([reflection; transmission], "reflection", "transmission")
title("Reflection and transmission coefficients")
xlabel("\theta_i [degree]")
xlim([0,90])
ylabel("Coefficient - module")


%% Energy conservation

theta = 0 : 90;
figure(2); 
grid on
hold on
S1z = plot(theta, abs(S1z_module),'LineWidth', 5);
S3z = plot(theta, abs(S3z_module),'LineWidth', 2, "Color", "magenta");
diff = plot(theta, abs(S1z_module) - abs(S3z_module) ,'LineWidth', 3, "Color", "black");
legend([S1z; S3z; diff], "|S1z|", "|S3z|", "|S1z|+|S3z|")
title("Conservation of energy")
xlabel("\theta_i [degree]")
xlim([0,90])
ylabel("|Sz|")


