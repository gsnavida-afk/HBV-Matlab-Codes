function Y = model_fit(x,year)
beta = x(1);
r1 = x(2);



b = [2.4662e-02,2.204e-02]; % estimated
mu=[1/43.26650546, 1/41.6608984]; % estimated
d=0.001; % estimated
epsilon=0.16;
k=4;
r2=0.025;
r3 = 0.01;
p=0.05;
q=0.1;
v=0.03;
w=0.05;

I_data1 = [6188 2874 2628 1982 2242 1745 1548]';
C_data1 = [424359 439848 453686 436344 444763 440314 445513]';


N0 = 51230704;
I0 = 6188;
L0 = I0/4;
C0 = 424359;
R0 = 0.3*N0;
S0 = N0-L0-I0-C0-R0;

init = [S0; L0; I0; C0; R0];
tspan = (2017:1:2023);

options = [];
[~,HB_model] = ode45(@model,year,init,options,b,mu,d,beta,epsilon,k,r1,r2,r3,w,p,q,v);
Y1 = HB_model(:,3);
Y2 = HB_model(:,4);
Y = [Y1 Y2];
