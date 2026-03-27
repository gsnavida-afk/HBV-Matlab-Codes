change=1:2:14;

b=2.7207e-02;
mu=1/43.42676617;
d=0.01;
beta=0.6503;
epsilon=0.16;
k=4;
r1=1.5381;
r2=0.025;
r3=0.01;
w=0.05;
p=0.05;
q=0.1;
v=0.11;
dummy=0.5;


L0 = 6446;
I0 = 25784;
C0 = 353403;
R0 = 16626604;
N0 = 49879811.5;
S0 = N0-L0-I0-C0-R0;


PRCC_var={'b','mu','d','beta','epsilon','k','r1','r2','r3','w','p','q','v','dummy'};

t_end=14;
tspan=(1:1:14);
time_points=1:2:t_end;

y0=[S0; L0; I0; C0; R0];

y_var_label={'S','L','I','C','R'};
