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
E2_r = sym ("E2_r");
M2_r = 5; 
h = sym("h");  % from cm to m
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

end

%% Choose two values for the incident angle theta 15 degrees & 45 degrees

gamma_1_true_value = 0.591031;
gamma_2_true_value = 0.732039;

gamma_1 = abs(Bs(16,1)); % theta+1 - Matlab indexing
gamma_2 = abs(Bs(46,1)); % theta+1 - Matlab indexing

%% Search space

E2r_values = linspace(0.8, 1.2, 21); % 1
H_values = linspace(0.008, 0.012, 21); % 0.01

for i = 1:length(E2r_values)
    for j = 1:length(H_values)
        E2r = E2r_values(i);
        H = H_values(j);
        expr = (gamma_1 - gamma_1_true_value) + (gamma_2 - gamma_2_true_value);
        substituted_expr = subs(expr, [E2_r, h], [E2r, H]);
        cost(i, j) = 1 / substituted_expr;
        
    end
end

%% Plot 

[E2rGrid, HGrid] = meshgrid(E2r_values, H_values);  

figure;  
mesh(E2rGrid, HGrid, abs(cost));  
xlabel('E2r');  
ylabel('H');  
zlabel('Cost');  
title('Mesh Plot of Cost Function'); 




