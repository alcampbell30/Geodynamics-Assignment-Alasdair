%***** RUN 2D MODEL FROM IMAGE ***********************************
% clear workspace
clear all; close all; %clc;

%Declare whether or not this is a test run. Switching this on will switch
%to constant coefficient variant. (YES/NO)
TEST = 'NO';
%Declare whether or not we wish to show that the total thermal energy in the system
%is constant in the test case (YES/NO)
conservation = 'NO';

% load model setup from image, interpolate to target grid size
W = 16e3; % domain width (must correspond to width of image) [m]
Nx = 400; % target grid size z-direction (select 100,200,400)
h = W/Nx; % grid spacing based on image width and target grid size
n_units = 9; % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);


% calculate sediment parameters
% set variables (values found online)

rho_particle = 2650; % density of sediments
rho_air = 1.225; % density of air
sa_phi = 0.23; % sand porosity proportion
gr_phi = 0.34; % gravel porosity proportion
si_phi = 0.15; % silt porosity proportion

% calculate bulk density of the sediment mass and air in the pores
rho_sand = (1-sa_phi) * rho_particle + (sa_phi * rho_air);
rho_grav = (1-gr_phi) * rho_particle + (gr_phi * rho_air);
rho_silt = (1-si_phi) * rho_particle + (si_phi * rho_air);


% material properties for each rock unit (update based on your calibration)
matprop = [
% unit conductivity(kT) density(rho0)   sp. heat capacity(Cp)   heat production(Hr)
1       3.678               2697.6        1000                  4.172e-6              %HE1
2       2.467               2700          874.5                 2.9e-6                %Gneiss 
3       3.218               2703.5        1000                  5.575e-6              %HE2
4       0.25                rho_sand      932                   1e-6                  %Sand
5       0.65                rho_grav      566                   1e-6                  %Gravel 
6       1.3                 2091.8        878                   1e-6                  %Clay (sea)
7       0.39                rho_silt      1088                  1e-6                  %Silt
8       0.61                1860          1510                  1e-6                  %Mud (Sea) (average 2512 Cp for wet mud and values for silt and sand)
9       0.024               rho_air       1000                  0];                   %air/water

  
% get coefficient fields based on spatial distribution of rock units from image
switch TEST
    case 'YES'
        rho0 = 2400*ones(Nz,Nx);
        Cp = 1000*ones(Nz,Nx);
        sigma = 1000*ones(Nz,Nx);
        Hr = ones(Nz, Nx); 
        air = units == 9;%make a variable to represent the air: i.e. no heat production, constant temp
        t_gradient = 35/1000;%Set a geothermal gradient of +30 degrees per 1000m of depth
        geotherm = [0, t_gradient]; 
        
        % pay attention if any unit conversion is required!
    case 'NO'
       rho0 = reshape(matprop(units,3),Nz,Nx); %isolate the third column in the table of data to isolate rho (Density)
       Cp = reshape(matprop(units,4),Nz,Nx);  %isolate the fourth column in the table of data to isolate Cp (Heat Capacity)
       sigma = reshape(matprop(units,2),Nz,Nx);  %isolate the second column in the table of data to isolate kT (Conductivity)
       Hr = reshape(matprop(units,5),Nz,Nx);  %isolate the fifth column in the table of data to isolate Hr (Radiogenic Heating)
       air = units == 9;%make a variable to represent the air: i.e. no heat production, constant temp
       t_gradient = 35/1000;%Set a geothermal gradient of +30 degrees per 1000m of depth
       geotherm = [0, t_gradient]; 
end

a = rho0.*Cp;
kT0 = sigma*1000./a; 

% continue setting remaining model parameters, then call model routine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ttop = 8; % surface temperature
Tbot = Ttop + D*geotherm(2); %Temperature at the bottom of the domain using geotherm gradient
cT = 1e-9; % kT T-dependence prefactor
mT = 2; % kT T-dependence powerlaw
g0 = 9.8; % gravity [m/s2]
aT = 1e-4; % thermal expansivity [1/C]

yr = 60*60*24*365.25; % seconds per year [s]
tend = 1e6*yr; % stopping time [s] 
CFL = 0.9; % Time step limiter 
nop = 100; % output figure produced every 'nop' steps

dTdto = 0;

switch conservation
    case 'YES' 

        run('./temperature_conservation_test.m');

    case 'NO'

         run('./helmsRK4.m');
end 




