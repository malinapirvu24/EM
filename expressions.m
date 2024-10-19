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
E1 = sym('E1'); 
k1 = sym('k1');

% second medium - layer
E2_r = 1;
M2_r = 5; 
h = sym('h');  
E2 = sym('E2');
k2 = sym('k2');

% third medium - substrate
E3_r = 1.25;
M3_r = 1; 
E3 = sym('E3');
k3 = sym('k3');

% define the ratio between first and third medium
ratio = sym('ratio');


%% Compute solution

theta = sym('theta');

k1z = k1*cos(theta);
k2z = sqrt(k2^2 - (k1 * sin(theta))^2);
k3z = sqrt(k3^2 - (k1 * sin(theta))^2);

% solution in Matrix form
A  = [1,          0,                       -1,                 -1;  ...
    0,      -exp(-1i*k3z*h),       exp(1i*k2z*h),         exp(- 1i*k2z*h); ...
    k1z/E1,       0,                      -k2z/E2,              k2z/E2; ...
    0,     (k3z/E3)*exp(-1i*k3z*h),  (k2z/E2)*exp(1i*k2z*h), (-k2z/E2)*exp(-1i*k2z*h)];

C = [-1; 0; k1z/E1; 0];

% compute solution for each incident angle -> assume Hi=1
A_inv = inv(A);
B = A_inv * C;


% get expression of reflection and transmission coefficients
gamma = B(1,1);
T = B(2,1)*ratio;





