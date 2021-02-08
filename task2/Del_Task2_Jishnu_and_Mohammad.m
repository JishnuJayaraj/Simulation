clear
close all
clc

% Given
h=0.0005; % delta x=delta y= h
Re=10000; % Reynolds number

x=0:h:1; % x axis
y=0:h:2*5/Re^0.5; % y axis

% Boundary Conditions
v=zeros(length(y),length(x)); % vertical velocity domain initialization
u=v; u(length(y),:)=1; u(:,1)=1; u(1,1)=0; % horizontal velocity domain initialization

A=zeros(length(y)-2,length(y)); % coefficient matrix initialization
r=zeros(length(y)-2,1); % result vector initialization

% Solving for u and v
for i=2:length(x)
    for j=2:length(y)-1
        A(j-1,j)=1+2/Re/h/u(j,i-1); % filling values in the diagonal of the coefficient matrix
        A(j-1,j+1)=v(j,i-1)/2/u(j,i-1)-1/Re/h/u(j,i-1); % filling values just right to the diagonal of the coefficient matrix
        A(j-1,j-1)=-v(j,i-1)/2/u(j,i-1)-1/Re/h/u(j,i-1); % filling values just left to the diagonal of the coefficient matrix

        r(j-1)=u(j,i-1); % filling results vector
    end
    
    r(1)=r(1)-A(1,1)*u(1,i); % correcting result vector values for lower boundary nodes
    r(end)=r(end)-A(end,end)*u(end,i); % correcting result vector values for upper boundary nodes
    M=A(:,2:end-1); % creating a new matrix as a correction of the coefficient matrix for lower and upper boundary nodes
    
    u(2:length(y)-1,i)=M\r; % solving the system of equations to get the horizontal velocity values for one level of nodes (one column) to be used in calculating the current level of vertical velocities and the next level of horizontal velocities
   
    for j=2:length(y)
        v(j,i)=v(j-1,i)-0.5*(u(j,i)-u(j,i-1)+u(j-1,i)-u(j-1,i-1)); % calculating vertical velocity values for one level of nodes (one column) to be used in calculating the next levels of horizontal and vertical velocities
    end 
end

% Showing results
SS=get(0,'ScreenSize');SW=SS(3);SH=SS(4); 
figure('NumberTitle','off','OuterPosition',[0 30 SW SH-30]);
subplot(2,1,1);imagesc(x,y,flipud(u)); colormap jet; colorbar; set(gca, 'XTick', linspace(0,max(x),21));set(gca, 'YTick', linspace(0,max(y),11),'YTickLabel', linspace(max(y),0,11)); grid; title('Distribution of horizontal dimentionless velocity on the whole domain'); xlabel('x'); ylabel('y'); 
subplot(2,1,2);imagesc(x,y,flipud(v)); colormap jet; colorbar; set(gca, 'XTick', linspace(0,max(x),21));set(gca, 'YTick', linspace(0,max(y),11),'YTickLabel', linspace(max(y),0,11)); grid; title('Distribution of vertical dimentionless velocity on the whole domain'); xlabel('x'); ylabel('y');

figure('NumberTitle','off','OuterPosition',[0 30 SW SH-30]);
subplot(2,2,1); plot(y,u(:,x==0.0005)*10); grid; title('Horizontal velocity porfile at x=0.0005m and U=10m/s'); xlabel('y'); ylabel('u (m/s)');
subplot(2,2,2); plot(y,v(:,x==0.0005)*10); grid; title('Vertical velocity profile at x=0.0005m and U=10m/s');  xlabel('y'); ylabel('v (m/s)');
subplot(2,2,3); plot(y,u(:,x==0.5)*10); grid; title('Horizontal velocity profile at x=0.5m and U=10m/s');  xlabel('y'); ylabel('u (m/s)');
subplot(2,2,4); plot(y,v(:,x==0.5)*10); grid; title('Vertical velocity profile at x=0.5m and U=10m/s');  xlabel('y'); ylabel('v (m/s)');

% Blasius exact solution 
deta=0.015;
eta=0:deta:8;
f=zeros(1,length(eta)); df=f; ddf=f; ddf(1)=0.3319;

for i=2:length(eta);
    f(i)=f(i-1)+df(i-1)*deta;
    df(i)=df(i-1)+ddf(i-1)*deta;
    ddf(i)=ddf(i-1)-0.5*f(i-1)*ddf(i-1)*deta;
end
uExact=df; 
vExact1=(1/Re/0.0005/2)^0.5*(eta.*df-f);
vExact2=(1/Re/0.5/2)^0.5*(eta.*df-f);

% Comparing numerical to exact solutions
figure('NumberTitle','off','OuterPosition',[0 30 SW SH-30]);
subplot(2,2,1); plot(y*(Re/0.0005)^0.5,u(:,x==0.0005),eta,uExact,'r'); xlim([0 max(eta)]);grid; title('The comparison between the exact and numerical horizontal velocity porfile at x=0.0005m'); xlabel('\eta = y/\delta(x)'); ylabel('u'); legend('Numerical Solution','Exact Solution','Location','SouthEast');
subplot(2,2,2); plot(y*(Re/0.0005)^0.5,v(:,x==0.0005),eta,vExact1,'r'); xlim([0 max(eta)]);grid; title('The comparison between the exact and numerical vertical velocity porfile at x=0.0005m'); xlabel('\eta  = y/\delta(x)'); ylabel('v'); legend('Numerical Solution','Exact Solution','Location','SouthEast');
subplot(2,2,3); plot(y*(Re/0.5)^0.5,u(:,x==0.5),eta,uExact,'r'); xlim([0 max(eta)]);grid; title('The comparison between the exact and numerical horizontal velocity porfile at x=0.5m'); xlabel('\eta  = y/\delta(x)'); ylabel('u'); legend('Numerical Solution','Exact Solution','Location','SouthEast');
subplot(2,2,4); plot(y*(Re/0.5)^0.5,v(:,x==0.5),eta,vExact2,'r');xlim([0 max(eta)]);grid; title('The comparison between the exact and numerical vertical velocity porfile at x=0.5m'); xlabel('\eta  = y/\delta(x)'); ylabel('v'); legend('Numerical Solution','Exact Solution','Location','SouthEast');
