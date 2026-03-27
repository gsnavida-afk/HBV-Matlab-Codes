clear all
clc


syms b mu d epsilon k r1 r2 r3 w p q v beta;

R01 = beta*b*k*w/((k+mu)*(r1+mu)*(r3+mu));
R02 = q + (epsilon*(p*q*r1+(r1+mu)*(1-q))/(d+r2+mu-(b*w*v)));
R0 = R01*R02;

x=[b mu d beta epsilon k r1 r2 r3 w p q v];
% simplify(diff(R0,b)*b/R0)
for i=1:13
    Sens(i)=simplify(diff(R0,x(i))*(x(i)/R0));
end
%%

b = 0.027207; % estimated
mu=0.02311; % estimated
d=0.01; 
beta=0.6503;
epsilon=20;
k=4;
r1=1.5381;
r2=0.025;
r3 = 0.01;
p=0.05;
q=0.1;
v=0.11;
w=0.05;


A=subs(Sens);
results=eval(A);


bar(results');
xlabel('Korean Model Parameters','FontSize',16); ylabel('Elasticity Index (R_{0}^{K})','FontSize',16);
%ylim([-1 1])
set(gca,'FontName','Times','XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13],...
    'XTickLabel',{'b','\mu','d','\beta','\epsilon','k','r_1','r_2','\sigma','w','p','q','v'} ...
    ,'FontSize',12,'YGrid','on','YMinorGrid','on')