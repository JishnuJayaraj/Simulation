function [ Ainv ] = Matrix_A(space,time)
No_of_Grid = 40; %%%%%%%%%Nuber of grid points
N_dim = No_of_Grid^2;
lmd = (time/(2*space^2));

A = zeros(N_dim,N_dim);
for i = 1:No_of_Grid
A(i,i) = 1;
end
for i = N_dim-No_of_Grid+1:N_dim
A(i,i) = 1;
end
A(No_of_Grid+1,No_of_Grid+1) = (1+4*lmd);
A(No_of_Grid+1,No_of_Grid+2) = -lmd;
A(No_of_Grid+1,(2*No_of_Grid)+1) = -lmd;
A(N_dim-No_of_Grid,N_dim-(2*No_of_Grid)) = -lmd;
A(N_dim-No_of_Grid,N_dim-No_of_Grid) = (1+4*lmd);
A(N_dim-No_of_Grid,N_dim-No_of_Grid-1) = -lmd;
for i = No_of_Grid+2:N_dim-No_of_Grid-1
A(i,i-1) = -lmd;
A(i,i) = (1+4*lmd);
A(i,i+1) = -lmd;
if(No_of_Grid+1<i && i<=(N_dim-No_of_Grid))
A(i,No_of_Grid+i) = -lmd;
end
if(No_of_Grid<i && i<=(N_dim-No_of_Grid))
A(i,i-No_of_Grid) = -lmd;
end
end
for i = 1:No_of_Grid
A(No_of_Grid*i,:)=0;
A(No_of_Grid*i,No_of_Grid*i)=1;
end
for i = 0:No_of_Grid-1
A(No_of_Grid*i+1,:)=0;
A(No_of_Grid*i+1,No_of_Grid*i+1)=1;
end
Ainv = A^-1;

end

