%***** RUN 2D MODEL FROM IMAGE ***********************************
% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size
W = 16e3; % domain width (must correspond to width of image) [m]
Nx = 200; % target grid size z-direction
h = W/Nx; % grid spacing based on image width and target grid size
n_units = 9; % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit (update based on your calibration)
matprop = [
% unit conductivity(kT) density(rho) heat capacity(Cp) heat production(Hr)
1       3.678               2697.6      1000                4.172   %HE1
2       1                   2000        1000                1       %Gneiss
3       1                   2000        1000                1       %Sand
4       3.218               2703.5      1000                5.575   %HE2
5       1                   2000        1000                1       %Gravel
6       1                   2000        1000                1       %Clay
7       1                   2000        1000                1       %Silt
8       1                   2000        1000                1       %Mud
9       1e-6                1000        1000                0];     % air/water
  
% get coefficient fields based on spatial distribution of rock units from image
% pay attention if any unit conversion is required!
rho = reshape(matprop(units,3),Nz,Nx); %isolate the third column in the table of data to isolate rho (Density)
Cp = reshape(matprop(units,4),Nz,Nx);  %isolate the fourth column in the table of data to isolate Cp (Heat Capacity)
kT = reshape(matprop(units,2),Nz,Nx);  %isolate the second column in the table of data to isolate kT (Conductivity)
Hr = reshape(matprop(units,5),Nz,Nx);  %isolate the fifth column in the table of data to isolate Hr (Radiogenic Heating)
air = units == 9;%make a variable to represent the air: i.e. no heat production, constant temp
geotherm = 30/1000; %Set a geothermal gradient of +30 degrees per 1000m of depth

%Test in constant coefficient unit value case (repeat for other variables)
%rho = 2400*ones(Nz,Nx);
% continue setting remaining model parameters, then call model routine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ttop = 5; % surface temperature
Tbot = Ttop + D*geotherm; %Temperature at the bottom of the domain using geotherm gradient
rho0 = 1000; % reference density [kg/m3]
kT0 = 1e-6; % Set an initial value of heat diffusivity that is to be used [m2/s]
cT = 1e-9; % kT T-dependence prefactor
mT = 2; % kT T-dependence powerlaw
g0 = 9.8; % gravity [m/s2]
aT = 1e-4; % thermal expansivity [1/C]

yr = 60*60*24*365; % seconds per year [s]
tend = 3e6*yr; % stopping time [s]
CFL = 3/4; % Time step limiter
nop = 100; % output figure produced every 'nop' steps
alpha = 0.99; % iterative step size limiter
beta = 0.95; % iterative lag parameter
tol = 1e-8; % residual tolerance
nup = 5; % update T, check residual every nup iterations
dTdto = 0;

run('./helmsdale_take_2.m');