function Y = model_fit2(x,year)
beta = x(1);
r = x(2);

b = 3.0324e-02; % estimated
mu=0.014318;
d=0.01;
k=4;
sigma=0.001;
v=0.11;
w=0.4263;


N0 = 106700000;
I0 = 40754;
L0 = I0/4;
R0 = 0.3*N0;
I02 = 40754;
S0 = N0-L0-I0-R0;

init = [S0; L0; I0; R0; I02];

options = [];
[~,HB_model] = ode45(@(t,x) ...
    modelaux(t,x,b,mu,d,beta,k,r,sigma,w,v),year,init,options);
Incidence = [HB_model(1,5); diff(HB_model(:,5))];

Y = Incidence;
end