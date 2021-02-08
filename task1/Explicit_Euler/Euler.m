%% EULER METHOD

% Solution of Heat equation in 2d using Explicit method

%% Parameters to define heat equation and define range in space and time

% x is length, y is width , t is time
% square grid x = y =(0 to 1)
L = 1; % length of x and y coordinates

%% Parameters to solve equation by Explicit method

% Since Explicit methods are stable only in small magnitude of dt, 
% given by the relation dt/(dx)^2, in the given set of time intervals
% only t=0.0001 is stable

maxt = 0.16;
dt = 0.0001 ;  % stable only when t = 0.0001
nt = maxt/dt;  % number of time step
nx = 41;       % no of grid points in x direction
ny = 41;       % no of grid points in y direction
dx = L/(nx-1);
dy = L/(ny-1);
b = dt/(dx)^2 ;

% to store the index of grid point at x = y =0.4
gridx = 0.4/dx ;  
gridy = 0.4/dy ;

%% Initial and Boundary conditions

t = 0:dt:0.16;
x = linspace(0,1,41);
y = linspace(0,1,41);


w = zeros(nx,ny);       % matrix to store all grid point values  
wnew = zeros(nx,ny);    % matrix to store grid point values for next iteration in loop

w(1,:)   = 1;
w(nx,:)  = 0;
w(:,1)   = 1-y;
w(:,nx)  = 1-y.*y;

wnew = w;

wpoint = zeros(nt,1);

%% Implimentation of explicit method
       % For each node points in n+1 time step, temperature is
       % calculated from previous time steps and stored in the matrix wnew
       
for time=0 : dt : 0.16
     
    for i=2:nx-1
       
        for  j=2:ny-1

         wnew(i,j)=w(i,j) + b*(w(i+1,j)+w(i-1,j)+w(i,j+1)+w(i,j-1) - 4* w(i,j));

         if i== gridx && j== gridy           % to get the values of point x=y=0.4
        
             wpoint(fix((time*10000)+1),1)=w(i,j);
             
         end
        end
    end    
    w = wnew;
    
    if time==0.01       % To plot temperature profile at time = 0.01
        figure(1);
        contourf(x,y,wnew);colorbar;title('2D Diffusion Heat equation')
        xlabel(' \it \bf x','FontSize',16)
        ylabel('\it \bf y','FontSize',16)
        title('Numerical solution to Heat equation at t=0.01','FontSize',15,'FontWeight','bold')
    end
    
    if time==0.02      % To plot temperature profile at time = 0.02
        figure(2);
        contourf(x,y,wnew);colorbar;
        xlabel(' \it \bf x','FontSize',16)
        ylabel('\it \bf y','FontSize',16)
        title('Numerical solution to Heat equation at t=0.02','FontSize',15,'FontWeight','bold')
    end
    
    if time==0.04     % To plot temperature profile at time = 0.04
        figure(3);
        contourf(x,y,wnew);colorbar;
        xlabel(' \it \bf x','FontSize',16)
        ylabel('\it \bf y','FontSize',16)
        title('Numerical solution to Heat equation at t=0.04','FontSize',15,'FontWeight','bold')
    end
    
    if time==0.08     % To plot temperature profile at time = 0.08
        figure(4);
        contourf(x,y,wnew);colorbar;
        xlabel(' \it \bf x','FontSize',16)
        ylabel('\it \bf y','FontSize',16)
        title('Numerical solution to Heat equation at t=0.08','FontSize',15,'FontWeight','bold')
    end
    
    if time==0.16     % To plot temperature profile at time =0.16
            figure(5);
            contourf(x,y,wnew);colorbar;
            xlabel(' \it \bf x','FontSize',16)
            ylabel('\it \bf y','FontSize',16)
            title('Numerical solution to Heat equation at t=0.16','FontSize',15,'FontWeight','bold')
    end  
    
end

% Time evolution of temperature at x = y = 0.4%wpoint = zeros(nx,ny);  % matrix to store temperature values at x=y=0.4

figure(6);
plot(t',wpoint,'linewidth',2);
xlabel( 'Time')
ylabel('Temperature')
title('Time evolution at point x = y= 0.4')

% Vertical Temperature profile at x=0.4 and t =0.16

figure(7);
plot(y,w(:,17));
xlabel(' Time')
ylabel('Temperature')
title('Vertical temperature profile at x=0.4 and t=0.16')
