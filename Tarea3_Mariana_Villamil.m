% Homework 3 for the Data Analysis Course
% Week 8: Finite Differences 2

% Student name:
stname = 'Mariana Carolina Villamil Sastre';
% Student codigo:
codigo = '201413371';
disp(['This is the work of ' stname ' with codigo ' codigo])
disp('')
disp('Homework 3: Finite Differences 2D')

%%
%Pore pressure diffusion in 2D equation with an explicit finite difference scheme
clear
%Physical parameters
L       =   200;        %   Width   [m]     
H       =   200;        %   Height   [m] 
Diff    =   0.01;       %   Difussivity [m^2/s]
Ppress   =   1;          %   Background pore pressure        [MPa]
Spress  =   10;         %   Square pression                [MPa] 


% Numerical parameters
nx      =   200;          %   # gridpoints in x-direction
nz      =   200;          %   # gridpoints in z-direction
dx      =   L/(nx-1);     %   Spacing of grid in x-direction
dz      =   H/(nz-1);     %   Spacing of grid in z-direction
[x2d,z2d] = meshgrid(0:dx:L, 0:dz:H);  % create grid

% Compute stable timestep
dt         = min([dx,dz])^2/Diff/4;
day     =   3600*24;    %   # seconds per day
nt=day/dt 
% Setup initial linear pressure profile
P =  ones(size(z2d))*Ppress;
P(90:110,90:110) =  Spress;
time    =   0;

for n=1:nt

    % Compute new temperature 
    Pnew = zeros(nz,nx);
    sx  = Diff*dt/dx^2;
    sz  = Diff*dt/dz^2;
    for j=2:nx-1
        for i=2:nz-1
            Pnew(i,j) = P(i,j) + (sx*(P(i+1,j)-(2*P(i,j))+P(i-1,j))) + ...
                (sz*(P(i,j+1)-(2*P(i,j))+P(i,j-1)));
        end
    end
    % Set boundary conditions
%    Pnew(1,:)   =  P(1 ,: ); % Dirichlet on top and bottom of grid
%    Pnew(nz,:)  =  P(nz,: );
%    for i=2:nz-1 % Neumann on left and right of grid
%        Pnew(i,1) = P(i,1) + (sx*( (2*P(i,2)) - (2*(P(i,1) + 2*dx*0)) )) + ...
%            (sz*( P(i+1,1) - (2*P(i,1)) + P(i-1,1) ));
%        Pnew(i,nx) = P(i,nx) + (sx*( (2*P(i,nx-1)) - (2*(P(i,nx) + 2*dx*0)) )) + ...
%            (sz*( P(i+1,nx) - (2*P(i,nx)) + P(i-1,nx) ));
%    end
    
     %%%% For Neumann on all borders
     Pnew(1,1) = P(1,2) - (2*dx*0);
     Pnew(nz,nx) = P(nz,nx-1) - (2*dx*0);
     Pnew(1,nx) = P(1,nx-1) - (2*dz*0);
     Pnew(nz,1) = P(nz-1,1) - (2*dz*0);
     for i=2:nz-1
         Pnew(i,1) = P(i,1) + (sx*( (2*P(i,2)) - (2*(P(i,1) + 2*dx*0)) )) + ...
             (sz*( P(i+1,1) - (2*P(i,1)) + P(i-1,1) ));
         Pnew(i,nx) = P(i,nx) + (sx*( (2*P(i,nx-1)) - (2*(P(i,nx) + 2*dx*0)) )) + ...
             (sz*( P(i+1,nx) - (2*P(i,nx)) + P(i-1,nx) ));
     end
     for i=2:nx-1
         Pnew(1,i) = P(1,i) + (sx*( (2*P(2,i)) - (2*(P(1,i) + 2*dx*0)) )) + ...
             (sz*( P(1,i+1) - (2*P(1,i)) + P(1,i-1) ));
         Pnew(nz,i) = P(nz,i) + (sx*( (2*P(nz-1,i)) - (2*(P(nz,i) + 2*dx*0)) )) + ...
             (sz*( P(nz,i+1) - (2*P(nz,i)) + P(nz,i-1) ));
     end
    
    P           =   Pnew;
    time        =   time+dt;
    % Plot solution every 50 timesteps
    if (mod(n,50)==0)
        figure(1), clf
        pcolor(x2d,z2d,Pnew); shading interp, colorbar
        hold on
        contour(x2d,z2d,Pnew,'k');
        xlabel('x [km]')
        ylabel('z [km]')
        zlabel('Temperature [^oC]')
        title(['Pressure evolution after ',num2str(time),'s'])
        drawnow
    end
end
CLS = Diff*dt/dx^2;
disp('The stability condition is reached')

%% Pore diffusion for 1 hour with diffusivity variable.
clear
%Physical parameters
L       =   200;        %   Width   [m]     
H       =   200;        %   Height   [m] 
Diff    =   0.01;       %   Difussivity [m^2/s]
Ppress   =   1;          %   Background pore pressure        [MPa]
Spress  =   10;         %   Square pression                [MPa] 


% Numerical parameters
nx      =   200;          %   # gridpoints in x-direction
nz      =   200;          %   # gridpoints in z-direction
dx      =   L/(nx-1);     %   Spacing of grid in x-direction
dz      =   H/(nz-1);     %   Spacing of grid in z-direction
[x2d,z2d] = meshgrid(0:dx:L, 0:dz:H);  % create grid

% Compute stable timestep
dt         = min([dx,dz])^2/0.1/4; 
hour=60*60 %Seconds per hour
nt=hour/dt;

% Setup initial linear pressure profile
P       =   ones(size(z2d))*Ppress;

P(90:110,90:110) =  Spress;

%Difussivity matrix
Dfs=   ones(size(z2d))*Diff;
Dfs(80:120,:) =  0.1;
time    =   0;

for n=1:nt

    % Compute new pressure with diffusivity average in 4 nodes 
    Pnew = zeros(nz,nx);

    for j=2:nx-1
        for i=2:nz-1
           
            sx  = (Dfs(i-1,j-1)+Dfs(i+1,j+1)+Dfs(i+1,j-1)+Dfs(i-1,j+1))*0.25*dt/dx^2;
            sz  = (Dfs(i-1,j-1)+Dfs(i+1,j+1)+Dfs(i+1,j-1)+Dfs(i-1,j+1))*0.25*dt/dz^2;
            
            Pnew(i,j) = P(i,j) + (sx*(P(i+1,j)-(2*P(i,j))+P(i-1,j))) + ...
                (sz*(P(i,j+1)-(2*P(i,j))+P(i,j-1)));
        end
    end
    
    % Set boundary conditions
    Pnew(1,:)   =  P(1 ,: ); % Dirichlet on top and bottom of grid
    Pnew(nz,:)  =  P(nz,: );
    for i=2:nz-1 % Neumann on left and right of grid
            sx  = Dfs(i,1)*dt/dx^2;
            sz  = Dfs(i,1)*dt/dz^2;
        Pnew(i,1) = P(i,1) + (sx*( (2*P(i,2)) - (2*(P(i,1) + 2*dx*0)) )) + ...
            (sz*( P(i+1,1) - (2*P(i,1)) + P(i-1,1) ));
        Pnew(i,nx) = P(i,nx) + (sx*( (2*P(i,nx-1)) - (2*(P(i,nx) + 2*dx*0)) )) + ...
            (sz*( P(i+1,nx) - (2*P(i,nx)) + P(i-1,nx) ));
    end
    
    P           =   Pnew;
    P(90:110,90:110) =  Spress;
    time        =   time+dt;
    
    % Plot solution every 50 timesteps
    if (mod(n,50)==0)
        figure(1), clf
        pcolor(x2d,z2d,Pnew); shading interp, colorbar
        hold on
        contour(x2d,z2d,Pnew);
        xlabel('x [m]')
        ylabel('z [m]')
        zlabel('Pressure [^MPa]')
        title(['Pressure evolution after ',num2str(time),' s'])
        drawnow
    end
end
CLS = Diff*dt/dx^2;
%% Pore diffusion for 1 day with diffusivity variable.
clear
%Physical parameters
L       =   200;        %   Width   [m]     
H       =   200;        %   Height   [m] 
Diff    =   0.01;       %   Difussivity [m^2/s]
Ppress   =   1;          %   Background pore pressure        [MPa]
Spress  =   10;         %   Square pression                [MPa] 


% Numerical parameters
nx      =   200;          %   # gridpoints in x-direction
nz      =   200;          %   # gridpoints in z-direction
dx      =   L/(nx-1);     %   Spacing of grid in x-direction
dz      =   H/(nz-1);     %   Spacing of grid in z-direction
[x2d,z2d] = meshgrid(0:dx:L, 0:dz:H);  % create grid

% Compute stable timestep
dt         = min([dx,dz])^2/0.2/4; 
day=3600*24 %Seconds per day
nt=day/dt;

% Setup initial linear pressure profile
P       =   ones(size(z2d))*Ppress;

P(90:110,90:110) =  Spress;

%Difussivity matrix
Dfs=   ones(size(z2d))*Diff;
Dfs(80:120,:) =  0.1;
time    =   0;

for n=1:nt

    % Compute new pressure with diffusivity average in 4 nodes 
    Pnew = zeros(nz,nx);

    for j=2:nx-1
        for i=2:nz-1
           
            sx  = (Dfs(i-1,j-1)+Dfs(i+1,j+1)+Dfs(i+1,j-1)+Dfs(i-1,j+1))*0.25*dt/dx^2;
            sz  = (Dfs(i-1,j-1)+Dfs(i+1,j+1)+Dfs(i+1,j-1)+Dfs(i-1,j+1))*0.25*dt/dz^2;
            
            Pnew(i,j) = P(i,j) + (sx*(P(i+1,j)-(2*P(i,j))+P(i-1,j))) + ...
                (sz*(P(i,j+1)-(2*P(i,j))+P(i,j-1)));
        end
    end
    
    % Set boundary conditions
    Pnew(1,:)   =  P(1 ,: ); % Dirichlet on top and bottom of grid
    Pnew(nz,:)  =  P(nz,: );
    for i=2:nz-1 % Neumann on left and right of grid
            sx  = Dfs(i,1)*dt/dx^2;
            sz  = Dfs(i,1)*dt/dz^2;
        Pnew(i,1) = P(i,1) + (sx*( (2*P(i,2)) - (2*(P(i,1) + 2*dx*0)) )) + ...
            (sz*( P(i+1,1) - (2*P(i,1)) + P(i-1,1) ));
        Pnew(i,nx) = P(i,nx) + (sx*( (2*P(i,nx-1)) - (2*(P(i,nx) + 2*dx*0)) )) + ...
            (sz*( P(i+1,nx) - (2*P(i,nx)) + P(i-1,nx) ));
    end
    
    P           =   Pnew;
    P(90:110,90:110) =  Spress;
    time        =   time+dt;
    
    % Plot solution every 50 timesteps
    if (mod(n,1000)==0)
        figure(1), clf
        pcolor(x2d,z2d,Pnew); shading interp, colorbar
        hold on
        contour(x2d,z2d,Pnew);
        xlabel('x [m]')
        ylabel('z [m]')
        zlabel('Pressure [^MPa]')
        title(['Pressure evolution after ',num2str(time),' s'])
        drawnow
    end
end


CLS = Diff*dt/dx^2;
