clear all
clc

load grids     
m = [0.4,0.6,0.8,1.2,1,0.1];     % the innput vector (configuration)
[u,T] = FE(medium,m);
plot_solution(medium,u)          % uncomment and run to see plots