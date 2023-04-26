function [u,T] = FE(tri,conf)

k1 = conf(1);
k2 = conf(2);
k3 = conf(3);
k4 = conf(4);
k5 = conf(5);
Bi = conf(6);

ki = [k1,k2,k3,k4,k5];

n = tri.nodes;
x = zeros(3,3);
c = zeros(3,3);
Alocal = zeros(3,3);
A = sparse(n,n);
F = sparse(n,1);

for r = 1:5                      %loop through the regions
for k = 1:size(tri.theta{r},1)   %loop through elements
        v1 = tri.theta{r}(k,1);  
        v2 = tri.theta{r}(k,2);
        v3 = tri.theta{r}(k,3);
        
        x(1,1) = 1;
        x(2,1) = 1;
        x(3,1) = 1;
        x(1,2) = tri.coor(v1,1);
        x(2,2) = tri.coor(v2,1);
        x(3,2) = tri.coor(v3,1);
        x(1,3) = tri.coor(v1,2);
        x(2,3) = tri.coor(v2,2);
        x(3,3) = tri.coor(v3,2);
        c = inv(x);                %finding the coeffs for H
        dA = 0.5*abs(det(x));      %area of triangular element
      
        
%%%%%%%%%%%% Domain Interior %%%%%%%%%%%% 
    for a = 1:3  
        for b =1:3
           Alocal(a,b) = ki(r)*(c(2,a)*c(2,b)+c(3,a)*c(3,b))*dA;   
           i = tri.theta{r}(k,a);
           j = tri.theta{r}(k,b);
           A(i,j) = A(i,j) + Alocal(a,b);
        end
    end
end
end
    
%%%%%%%%%%%% Boundaries except for Root %%%%%%%%%%
for k = 1: size(tri.theta{6},1)
    n1 = tri.theta{6}(k,1);
    n2 = tri.theta{6}(k,2);
    x1 = tri.coor(n1,1);
    y1 = tri.coor(n1,2);
    x2 = tri.coor(n2,1);
    y2 = tri.coor(n2,2);
    h = ((x1-x2)^2 + (y1-y2)^2)^0.5;
    Alocal = Bi*(h/3)*[1 0.5; 0.5 1];
    A(n1,n1) = A(n1,n1) + Alocal(1,1);
    A(n1,n2) = A(n1,n2) + Alocal(1,2);
    A(n2,n1) = A(n2,n1) + Alocal(2,1);
    A(n2,n2) = A(n2,n2) + Alocal(2,2);
end 


%%%%%%%%%%%%%%%% Root Boundary %%%%%%%%%%%%%
for k = 1: size(tri.theta{7},1)
    n1 = tri.theta{7}(k,1);
    n2 = tri.theta{7}(k,2);
    x1 = tri.coor(n1,1);
    y1 = tri.coor(n1,2);
    x2 = tri.coor(n2,1);
    y2 = tri.coor(n2,2);
    h = ((x1-x2)^2 + (y1-y2)^2)^0.5;
    Flocal = (h/2)*[1;1];
    F(n1) = F(n1) + Flocal(1);
    F(n2) = F(n2) + Flocal(2);
end 

u  = A\F;    %%%% solve for u
T= u'*F;     %%%% solve for T

end


