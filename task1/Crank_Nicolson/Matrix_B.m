function [ Bn ] = Matrix_B(space,time )

No_of_Grid = 40; %%%%%%%%%Nuber of grid points
N_dim = No_of_Grid^2;
lmd = (time/(2*space^2));
B = zeros(N_dim,N_dim);
for i = 1:No_of_Grid
B(i,i) = 1;
end
for i = N_dim-No_of_Grid+1:N_dim
B(i,i) = 1;
end
B(No_of_Grid+1,No_of_Grid+1) = (1-4*lmd);
B(No_of_Grid+1,No_of_Grid+2) = lmd;
B(No_of_Grid+1,(2*No_of_Grid)+1) = lmd;
B(N_dim-No_of_Grid,N_dim-(2*No_of_Grid)) = lmd;
B(N_dim-No_of_Grid,N_dim-No_of_Grid) = (1-4*lmd);
B(N_dim-No_of_Grid,N_dim-No_of_Grid-1) = lmd;
for i = No_of_Grid+2:N_dim-No_of_Grid-1
B(i,i-1) = lmd;
B(i,i) = (1-4*lmd);
B(i,i+1) = lmd;
if(No_of_Grid+1<i && i<=(N_dim-No_of_Grid))
B(i,No_of_Grid+i) = lmd;
end
if(No_of_Grid<i && i<=(N_dim-No_of_Grid))
B(i,i-No_of_Grid) = lmd;
end
end
for i = 1:No_of_Grid
B(No_of_Grid*i,:)=0;
B(No_of_Grid*i,No_of_Grid*i)=1;
end
for i = 0:No_of_Grid-1
B(No_of_Grid*i+1,:)=0;
B(No_of_Grid*i+1,No_of_Grid*i+1)=1;
Bn=B;
end
