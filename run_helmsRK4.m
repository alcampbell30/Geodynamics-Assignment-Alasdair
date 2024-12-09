%***** RUN 2D MODEL FROM IMAGE ***********************************
% clear workspace
clear all; close all; %clc;

%Declare whether or not this is a test run. Switching this on will switch to constant coefficient variant. 
TEST = 'NO';
NN = [100,200,400];

for nn = 1:3

% load model setup from image, interpolate to target grid size
W = 16e3; % domain width (must correspond to width of image) [m]
Nx = NN(nn); % target grid size z-direction
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
% unit conductivity(kT) density(rho0) sp. heat capacity(Cp)   heat production(Hr)
1       3.678               2697.6      845                   4.172e-6              %HE1
2       2.467               2750        775                   5.7e-6                %Gneiss 
3       3.218               2703.5      845                   5.575e-6              %HE2
4       0.272               rho_sand    830                   1e-6                  %Sand
5       1.075               rho_grav    1000                  1e-6                  %Gravel
6       1.3                 2000        1381                  1e-6                  %Clay (sea)
7       2.49                rho_silt    1000                  1e-6                  %Silt
8       0.61                2000        2512                  1e-6                  %Mud (Sea)
9       1e-6                rho_air     1012                  0];                   %air/water

  
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
Ttop = 10; % surface temperature
Tbot = Ttop + D*geotherm(2); %Temperature at the bottom of the domain using geotherm gradient
%rho0 = 1000; % reference density [kg/m3]
%kT0 = 1e-6; % Set an initial value of heat diffusivity that is to be used [m2/s]
cT = 1e-9; % kT T-dependence prefactor
mT = 2; % kT T-dependence powerlaw
g0 = 9.8; % gravity [m/s2]
aT = 1e-4; % thermal expansivity [1/C]

yr = 60*60*24*365.25; % seconds per year [s]
tend = 1e6*yr; % stopping time [s]
CFL = 0.9; % Time step limiter 
nop = 100; % output figure produced every 'nop' steps

dTdto = 0;

run('./helmsRK4.m');
end

%Plot numerical error
Ex(nn)  = Errx;
Dh(nn) = h;

%plot numerical error in x 
figure(); 
loglog(Dh,Ex,'ro','LineWidth',1.5,'MarkerSize',8); axis tight; box on; hold on
loglog(Dh,Ex(1).*[1,1/2,1/4].^1,'k-','LineWidth',0.7)
loglog(Dh,Ex(1).*[1,1/2,1/4].^2,'k-','LineWidth',0.9)
loglog(Dh,Ex(1).*[1,1/2,1/4].^3,'k-','LineWidth',1.1)
loglog(Dh,Ex(1).*[1,1/2,1/4].^4,'k-','LineWidth',1.3)
loglog(Dh,Ex(1).*[1,1/2,1/4].^5,'k-','LineWidth',1.5)
xlabel('Step size','FontSize',18)
ylabel('Numerical error','FontSize',18)
title('Numerical Convergence in Space','FontSize',20)

Ez(nn)  = Errz;
Dt(nn) = t;

%plot numerical error in z
figure(); 
loglog(Dt,Ez,'ro','LineWidth',1.5,'MarkerSize',8); axis tight; box on; hold on
loglog(Dt,Ez(1).*[1,1/2,1/4].^1,'k-','LineWidth',0.7)
loglog(Dt,Ez(1).*[1,1/2,1/4].^2,'k-','LineWidth',0.9)
loglog(Dt,Ez(1).*[1,1/2,1/4].^3,'k-','LineWidth',1.1)
loglog(Dt,Ez(1).*[1,1/2,1/4].^4,'k-','LineWidth',1.3)
loglog(Dt,Ez(1).*[1,1/2,1/4].^5,'k-','LineWidth',1.5)
xlabel('Step size','FontSize',18)
ylabel('Numerical error','FontSize',18)
title('Numerical Convergence in Space','FontSize',20)
