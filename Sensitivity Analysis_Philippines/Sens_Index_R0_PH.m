clear all
clc


syms b mu d  k r sigma w p q v beta;

R0 = beta*b*k*w/((k+mu)*(sigma+mu)*(d+r+mu-b*w*v));

x=[b mu d beta k r sigma w v];
% simplify(diff(R0,b)*b/R0)
for i=1:9
    Sens(i)=simplify(diff(R0,x(i))*(x(i)/R0));
end
%%

b = 3.0324e-02;
mu=0.014318;
d=0.01;
beta=1.5045;
k=4;
r=0.8626;
sigma=0.001;
v=0.11;
w=0.4263;


A=subs(Sens);
results=eval(A);


bar(results');
set(gca,'FontName','Times','XTick',[1 2 3 4 5 6 7 8 9],...
    'XTickLabel',{'b','\mu','d','\beta','k','\rho','\sigma','w','v'} ...
    ,'FontSize',12,'YGrid','on','YMinorGrid','on')