%***** 1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT *******************
%***** Initialise Model Setup

%Image data obtained
% create x-coordinate vectors
xc = h/2:h:W-h/2; % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2; % z-coordinate vector for cell centre positions [m]
xf = 0:h:W; % x-coordinate vectore for cell face positions [m]
zf = 0:h:D; % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc); % create 2D coordinate arrays
%Coordinate grid on the output / for code interpretation: use this as the
%model grid

%h = W/Nx;
%xc = h/2:h:W-h/2;
%zc = h/2:h:D-h/2;
%[Xc,Zc] = meshgrid(xc,zc);

%Insulating Sides
ix3 = [ 1,1:Nx,Nx ];
ix5 = [ 1, 1,1:Nx,Nx,Nx ];
ix7 = [ 1, 1, 1,1:Nx,Nx,Nx,Nx];

%Set up Insulating Top and Bottom
iz3 = [ 1,1:Nz,Nz ];
iz5 = [ 1, 1,1:Nz,Nz,Nz ];
iz7 = [ 1, 1, 1,1:Nz,Nz,Nz,Nz];

% create smooth random perturbation field
rng(15);
dr = randn(Nz,Nx);
for ii = 1:10
dr = dr + (diff(dr(iz3,:),2,1) + diff(dr(:,ix3),2,2))/8;
end

% set initial condition for temperature at cell centres
T = Ttop + geotherm.*Zc; % initialise T array on linear gradient
T(air) = Ttop;
% initialise density and mobility
rho = rho0.*(1 - aT.*(T-Ttop));
kT = kT0.*ones(Nz,Nx);
% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,0,yr)
%***** Solve Model Equations
dt = CFL * min((h/2)^2/max(kT(:))); % initial time step [s]
t = 0; % initial time [s]
k = 0; % initial time step count
dTdt = 0;

% loop through time steps until stopping time reached
while t <= tend
    % increment time and step count
    t = t+dt;
    k = k+1;
    % print time step header
    fprintf(1,'\n\n***** step = %d; dt = %1.3e; time = %1.3e \n\n',k,dt,t)
    % store old temperature and rate
    dTdto = dTdt;
    To = T;
    resnorm = 1; % initialise residual norm
    it = 0; % initialise iteration count
    % loop through pseudo-transient iterations until convergence criterion reached
    while resnorm > tol | it<10
    % update temperature every 'nup' iterations
        if ~mod(it,nup) && k>=1
        % get rate of change

        
        dTdt = diffusion(T,kT,h,ix3,iz3,geotherm, Hr, rho, Cp);
        % get temperature residual
        res_T = (T - To)/dt - (dTdt + dTdto)/2;
        % set isothermal boundaries on top/bot
        res_T(1 ,:) = 0;
        res_T(end,:) = 0;
        % get solution update
        upd_T = - alpha*res_T*dt/2;
        % keep air isothermal
        upd_T(air) = 0;
        % update solution
        T = T + upd_T;
        end
    
    it = it+1; % increment iteration count
    % get residual norm and print convergence every 'nup' iterations
    if ~mod(it,nup)
    resnorm = norm(upd_T(:),2)./norm(T(:)+eps,2);
        if isnan(resnorm); error('!!! Solver failed with nan !!!'); end
        fprintf(1,' it = %d; res = %e \n',it,resnorm);
        end
end

% plot model progress every 'nop' time steps
if ~mod(k,nop)
makefig(xc,zc,T,t,yr);
end
end
%***** Utility Functions ************************************************
% Function to make output figure
function makefig(x,z,T,t,yr)

    % plot temperature in subplot 1
    imagesc(x,z,T); axis equal; colorbar; hold on
    contour(x,z,T,[100,150,200],'k');
    drawnow
    ylabel('Depth [m]','FontSize',15)
    title(['Temperature; time = ',num2str(t/yr),' yr'],'FontSize',18)
end

% Function to calculate diffusion rate
function [dTdt] = diffusion(f,k,h,ix,iz,geotherm, Hr, rho, Cp)
% calculate heat flux by diffusion
    kx = (k(:,ix(1:end-1)) + k(:,ix(2:end)))/2; %find the heat conductivity on the cell faces / edges in x direction
    kz = (k(iz(1:end-1),:) + k(iz(2:end),:))/2; %find the heat conductivity on the cell faces / edges in z directionx

    qx = - kx .* diff(f(:,ix), 1, 2)/h; %Calculate heat flux across the cell boundaries using kx in the x direction
    qz = - kz .* diff(f(iz,:), 1, 1)/h; %Calculate heat flux across the cell boundaries using kz in the z direction 

    % basalt boundary
    qz(end,:) = - kz(end,:) .* geotherm;

    % calculate flux balance for rate of change
    dTdt_diffusion = - (diff(qx,1,2)/h+diff(qz,1,1)/h);

    % add heat source term using Hr data from Matprop table
    heat_source = Hr./(rho.*Cp); % Source term due to heat production.

    % total rate of change of temperature
    dTdt = dTdt_diffusion + heat_source;
end