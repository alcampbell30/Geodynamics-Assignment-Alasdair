%***** RUN 2D MODEL FROM IMAGE ***********************************
% clear workspace
clear all; close all; %clc;

%Declare whether or not this is a test run. Switching this on will switch to constant coefficient variant. 
TEST = 'NO';

% load model setup from image, interpolate to target grid size
W = 16e3; % domain width (must correspond to width of image) [m]
Nx = 200; % target grid size z-direction
h = W/Nx; % grid spacing based on image width and target grid size
n_units = 9; % number of rock units contained in image
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);


% calculate sediment parameters
% set variables (values found online)

rho_particle = 2650; % density of sediments
rho_air = 1.225; % density of air
sa_phi = 0.253; % sand porosity proportion
gr_phi = 0.325; % gravel porosity proportion
si_phi = 0.17; % silt porosity proportion

sa_particle_kT = 8; % sand particle thermal conductivity
gr_particle_kT = 5; % gravel particle thermal conductivity
si_particle_kT = 3.5; % silt particle thermal conductivity
air_kT = 0.025; % air thermal conductivity

% calculate bulk density of the sediment mass and air in the pores
rho_sand = (1-sa_phi) * rho_particle + (sa_phi * rho_air);
rho_grav = (1-gr_phi) * rho_particle + (gr_phi * rho_air);
rho_silt = (1-si_phi) * rho_particle + (si_phi * rho_air);


% material properties for each rock unit (update based on your calibration)
matprop = [
% unit conductivity(kT) density(rho0) heat capacity(Cp) heat production(Hr)
1       3.678               2697.6      845                 4.172                %HE1
2       1                   2800        1000                1                    %Gneiss
3       1                   rho_sand    830                 (1-sa_phi)*1         %Sand
4       3.218               2703.5      845                 5.575                %HE2
5       1                   rho_grav    1000                (1-gr_phi)*0.6       %Gravel
6       1                   1760        1381                2.75                 %Clay (Wet)
7       1                   rho_silt    1000                (1-si_phi)*1.4       %Silt
8       1                   2000        2512                3.5                  %Mud
9       1e-6                rho_air     1000                0];                  %air/water
  
% get coefficient fields based on spatial distribution of rock units from image
switch TEST
    case 'YES'
        rho0 = 2400*ones(Nz,Nx);
        Cp = 1000*ones(Nz,Nx);
        sigma = ones(Nz,Nx);
        Hr = ones(Nz, Nx); 
        air = units == 9;%make a variable to represent the air: i.e. no heat production, constant temp
        geotherm = 35/1000; %Set a geothermal gradient of +30 degrees per 1000m of depth
        
        % pay attention if any unit conversion is required!
    case 'NO'
       rho0 = reshape(matprop(units,3),Nz,Nx); %isolate the third column in the table of data to isolate rho (Density)
       Cp = reshape(matprop(units,4),Nz,Nx);  %isolate the fourth column in the table of data to isolate Cp (Heat Capacity)
       sigma = reshape(matprop(units,2),Nz,Nx);  %isolate the second column in the table of data to isolate kT (Conductivity)
       Hr = reshape(matprop(units,5),Nz,Nx);  %isolate the fifth column in the table of data to isolate Hr (Radiogenic Heating)
       air = units == 9;%make a variable to represent the air: i.e. no heat production, constant temp
       geotherm = 35/1000; %Set a geothermal gradient of +30 degrees per 1000m of depth

end

kT0 = sigma./rho0.*Cp;

% continue setting remaining model parameters, then call model routine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ttop = 5; % surface temperature
Tbot = Ttop + D*geotherm; %Temperature at the bottom of the domain using geotherm gradient
%rho0 = 1000; % reference density [kg/m3]
%kT0 = 1e-6; % Set an initial value of heat diffusivity that is to be used [m2/s]
cT = 1e-9; % kT T-dependence prefactor
mT = 2; % kT T-dependence powerlaw
g0 = 9.8; % gravity [m/s2]
aT = 1e-4; % thermal expansivity [1/C]

yr = 60*60*24*365; % seconds per year [s]
tend = 3e6*yr; % stopping time [s]
CFL = 0.95; % Time step limiter decreasing this will increase the time step, on short scale less accurate on long scale quicker
nop = 100; % output figure produced every 'nop' steps


dTdto = 0;

run('./helmsRK4.m');
