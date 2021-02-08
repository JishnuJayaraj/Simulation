clc
clear
close all

SS=get(0,'ScreenSize');SW=SS(3);SH=SS(4);

% main parameters
h=0.05;
dt=0.00001;
Re=[0.1 1];
Vel=1;
xy= 0:h:1;

% velocity matrix
u = zeros(length(xy),length(xy)+1);
v = zeros(length(xy)+1,length(xy));

uu=zeros(2,length(xy)+1);
vv=zeros(2,length(xy));

F = zeros(length(xy),length(xy)+1);
G = zeros(length(xy)+1,length(xy));

% pressure matrix
P = zeros(length(xy)+1,length(xy)+1);

% result vector
nElements=length(xy)-1;
s=nElements*nElements;

B = zeros(s,1);

% boundary conditions
F(1,:)=0;
F(length(xy),:)=0;

G(:,1)=0;
G(:,length(xy))=0;

% preparing the diagonal vectors for the coefficient matrix
nn=1:s;

d00=-4/h^2*ones(1,s);
d00(nn/nElements<1 | nn/nElements>nElements-1 | mod(nn,nElements)==0 | mod(nn,nElements)==1)=-3/h^2;
d00(1)=-2/h^2;
d00(end)=-2/h^2;
d00(nElements)=-2/h^2;
d00(s-nElements+1)=-2/h^2;

d01=1/h^2*ones(1,s-1);
d01(mod(1:s-1,nElements)==0)=0;

d20=1/h^2*ones(1,s-nElements);

% coefficient matrix
P1=sparse(diag(d00,0)+diag(d01,-1)+diag(d01,1)+diag(d20,nElements)+diag(d20,-nElements));

% time marching
for k=1:2
    for t=0:length(0:dt:0.1)
        
        % calculating F
        for i= 2: length(xy)-1
            for j= 2:length(xy)
                F(i,j)= u(i,j) + dt*((u(i+1,j) -2*u(i,j)+u(i-1,j)+u(i,j+1)-2*u(i,j)+u(i,j-1))*1/(Re(k)*h^2)-1/h*((u(i,j) +u(i+1,j)).^2/4-(u(i,j)+u(i-1,j)).^2/4+ (u(i,j+1)+u(i,j))*(v(i,j)+v(i+1,j))/4- (u(i,j)+u(i,j-1))*(v(i,j-1)+v(i+1,j-1))/4));
            end
        end
        
        % calulating G
        for i= 2: length(xy)
            for j= 2:length(xy)-1
                G(i,j)= v(i,j) + dt*((v(i+1,j) -2*v(i,j)+v(i-1,j)+v(i,j+1)-2*v(i,j)+v(i,j-1))*1/(Re(k)*h^2) -1/h*((v(i,j)+v(i,j+1)).^2/4-(v(i,j)+v(i,j-1)).^2/4 + (u(i,j)+u(i,j+1))*(v(i+1,j)+v(i,j))/4 - (u(i-1,j+1)+ u(i-1,j))*(v(i,j)+v(i-1,j))/4));
            end
        end
        
        index=1;
        for i=2:length(xy)
            for j=2: length(xy)
                B(index,1) = (1/(h*dt))*(F(i,j)-F(i-1,j)+G(i,j)-G(i,j-1));
                index=index+1;
            end
        end
        B(mod(1:s-1,nElements)==0 | mod(1:s-1,nElements)==1)=B(mod(1:s-1,nElements)==0 | mod(1:s-1,nElements)==1)-1/h^2; B(1)=0;
        
        P_result= P1\B;
        
        P(2:length(xy),2:length(xy))=reshape(P_result,length(xy)-1,length(xy)-1)';
        
        P(2:length(xy),1)= P(2:length(xy),2);
        P(2:length(xy),length(xy)+1)= P(2:length(xy),length(xy));
        P(1,2:length(xy))= P(2,2:length(xy));
        P(length(xy)+1,2:length(xy))= P(length(xy),2:length(xy));
        
        for i=2:length(xy)-1
            for j=2:length(xy)
                u(i,j) = F(i,j)-dt/h*(P(i+1,j)-P(i,j));
            end
        end
        
        for i=2:length(xy)
            for j=2:length(xy)-1
                v(i,j)= G(i,j)-dt/h*(P(i,j+1)-P(i,j));
            end
        end
        
        % updating boundary conditions for the new time step
        u(:,1)= -u(:,2);
        u(2:end-1,end)= 2*Vel-u(2:end-1,end-1); 
        
        v(1,:)= -v(2,:);
        v(end,:)= -v(end-1,:);
        
        % plotting resutls
        if((t==2001 || t== 4001 || t==6001 || t==8001 || t==length(0:dt:0.1)) && k==1)
            
            figure('Name',['Navier stokes equation for Pressure contour and Veloicity Streamline at t=',num2str((t-1)/100000)],'NumberTitle','off','OuterPosition',[0 30 SW SH-30])
            subplot(1,2,1);
            contourf(xy(2:end),xy(2:end),P(2:end-1,2:end-1)',50); colorbar
            xlabel(' Spatial co-ordinate (x)')
            ylabel('Spatial co-ordinate (y)')
            title(['Pressure contour plot for Navier stokes equation at t=',num2str((t-1)/100000)])
            
            subplot(1,2,2);
            quiver(xy(2:end),xy(2:end),u(2:end,2:end-1)',v(2:end-1,2:end)');
            streamline(xy(2:end),xy(2:end),u(2:end,2:end-1)',v(2:end-1,2:end)',h:0.2:1,h:0.2:1)
            axis([-0.1 1.1 -0.1 1.1])
            xlabel(' Spatial co-ordinate (x)')
            ylabel('Spatial co-ordinate (y)')
            title(['Streamline diagram for Navier stokes equation at t=',num2str((t-1)/100000)])
            
        end
        
    end
    
   %saving comparison vectors
   uu(k,:)=u(floor(length(xy)/2),:);
   vv(k,:)=v(floor((length(xy)+1)/2),:);
end

% calculating mean horizontal velocities to align them with the boundaries of the study domain
 uMean=uu;
 uMean(:,end)=[];
 uMean=[[0;0],uMean];
 uMean=(uMean+uu)./2;
 uMean(:,1)=[];

% plotting comparison
figure('Name','Horizontal and vertical velocity distributions at x=0.5 and t=0.1 for Re=0.1 and Re=1','NumberTitle','off','OuterPosition',[0 30 SW SH-30]);
subplot(1,2,1); 
plot(linspace(0,1,length(xy)),uMean(1,:)); hold on
plot(linspace(0,1,length(xy)),uMean(2,:));
legend('u @ Re=0.1','u @ Re=1','Location','NorthWest');
title('Comparing horizontal velocities at Re=0.1 and Re=1')
xlabel('Spatial co-ordinate (y)');
ylabel('Horizontal dimentionless velocity (u)');
grid

subplot(1,2,2); 
plot(xy,vv(1,:) ); hold on
plot(xy,vv(2,:) );
legend('v @ Re=0.1','v @ Re=1','Location','NorthWest');
title('Comparing vertical velocities at Re=0.1 and Re=1')
xlabel('Spatial co-ordinate (y)');
ylabel('Vertical dimentionless velocity (v)');
grid