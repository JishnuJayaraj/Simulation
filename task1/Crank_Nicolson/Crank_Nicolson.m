
% Crank-Nicolson method of solution

% Solving by Alternating Direction Implicit
%two-level solution method. This approach ensures that at each time step,
%there are no more than three unknowns to solve.


% Initial conditions

space = 0.025;
time= 0.0001;
T=0.16;
t=T/time;

% Parameters to solve equation by Crank Nicolson

No_of_Grid = 40;        % Number of grid points
N_dim = No_of_Grid^2;   
p = zeros(N_dim,1);     % Column matrix to be solved 
q=p;                    

% variable to get the co ordinates of point x=y=0.4

Data= zeros(t,1);
m=0.4/space;           
n=0.4/space;

% Boundary conditions initialisation

for i = 1:No_of_Grid-1         
    q(No_of_Grid*i,1) = 0;
    q(No_of_Grid*i+1,1) = 1;
   
end

for i = 1:No_of_Grid
    q(i,1) = 1-(i*space);
    q(N_dim-No_of_Grid+i,1) = (1-((i*space)^2));
end

inv = Matrix_A(space,time);     % Function to evaluate inverse of matrix

% Evolution with time

for i = 1:t
    p = inv * Matrix_B(space,time) * q ;
    q=p;
    Data(i,1) = p(m*n,1);
end

% Converting the solution to matrix form
c=zeros(No_of_Grid,No_of_Grid); 
for i = 1:No_of_Grid
    x = ((i-1)*No_of_Grid)+1;
    y = i*No_of_Grid;
    c(:,i)= p(x:y,1);
end

% Ploting the results

 figure(1)           % Ploting the contour of heat solution
 contourf(c,13);
 colorbar('Direction','reverse');
 xlabel('spatial coordinate x','FontWeight','bold');
 ylabel ('spatial coordinate y','FontWeight','bold');
 title('2D Diffusion Heat equation')

 figure(2)           % Plotting Temperature at point x=y=0.4 
 x = time:time:T;
 xlabel('Time Step','FontWeight','bold');
 ylabel ('Temperature at (0.4, 0.4)','FontWeight','bold');
 title('Time evolution of temperature at x=y=0.4 ');
 plot(x,Data,'--gs','MarkerFacecolor',[0.5 0.5 0.5],'MarkerEdgecolor','b','MarkerSize',2);
 grid on;
 hold on;

 figure(3)   % Plotting Temperature distribution on vertical line x=0.4 at t=0.16
 info1 = c(:,16);
 x1 = space:space:1;
 plot(x1,info1,'--gs','MarkerFacecolor',[0.5 0.5 0.5],'MarkerEdgecolor','b','MarkerSize',5);
 xlabel('y','FontWeight','bold');
 ylabel ('Temperature on line x= 0.4 at t=0.16','FontWeight','bold');
 title('Vertical Temperature Profile at x=0.4 and t=0.16 '); 
 grid on;
 hold on;

