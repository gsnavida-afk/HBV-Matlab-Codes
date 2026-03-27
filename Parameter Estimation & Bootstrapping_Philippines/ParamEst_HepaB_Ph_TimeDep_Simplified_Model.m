% Parameter Estimation

clc;clear; close all
load IncidenceHepaBPh1.mat %Hepa B prevalence data from 2010-2023

t_data = IncidenceHepaBPh1(:,1);
I_data = IncidenceHepaBPh1(:,2);
%% Optimization for the First Period (2017-2019)
init_par = [2.6797 0.025]; % initial values
LB = [0 0];
UB = [5 1];
options = optimoptions("lsqcurvefit",'StepTolerance',1e-15,'MaxFunctionEvaluations',1e+5);
xopt1= lsqcurvefit(@model_fit1,init_par,t_data(1:3),I_data(1:3),LB,UB,options);

beta = xopt1(1);
r = xopt1(2);

display(xopt1);
%% Best Fit Model for the First Period (2017-2019)

b = 3.0324e-02; % estimated
mu=0.014318;
d=0.01;
k=4;
sigma=0.001;
v=0.11;
w=0.4263;

N0 = 106700000;
I0 = 72598;
L0 = I0/4;
R0 = 0.3*N0;
I02 = 72598;
S0 = N0-L0-I0-R0;

init = [S0; L0; I0; R0; I02];

options = [];
[~,X] = ode45(@modelaux,t_data(1:3),init,options,b,mu,d,beta,k,r,sigma,w,v);
I1 = [X(1,5); diff(X(:,5))];

Rnot1 = beta*b*k*w/((k+mu)*(sigma+mu)*(d+r+mu-b*w*v));

display(Rnot1);



%% Optimization for the Second Period (2020-2023)
init_par = [2.6797 0.025]; % initial values
LB = [0 0];
UB = [5 1];
options = optimoptions("lsqcurvefit",'StepTolerance',1e-15,'MaxFunctionEvaluations',1e+5);
xopt2= lsqcurvefit(@model_fit2,init_par,t_data(4:end),I_data(4:end),LB,UB,options);

beta = xopt2(1);
r = xopt2(2);

display(xopt2);
%% Best Fit Model for the Second Period (2020-2023)

N0 = 106700000;
I0 = 40754;
L0 = I0/4;
R0 = 0.3*N0;
I02 = 40754;
S0 = N0-L0-I0-R0;

init = [S0; L0; I0; R0; I02];

options = [];
[~,X2] = ode45(@modelaux,t_data(4:end),init,options,b,mu,d,beta,k,r,sigma,w,v);
I2 = [X2(1,5); diff(X2(:,5))];

Rnot2 = beta*b*k*w/((k+mu)*(sigma+mu)*(d+r+mu-b*w*v));

display(Rnot2);

N = 114891199;
Snot = b*w*N/(sigma+mu);
Pnot = b*(mu*(1-w)+sigma)*N/(mu*(sigma+mu));

display(num2str([Snot 0 0 Pnot]));

T = b*w*v + ((k+mu)*(d+r+mu-b*w*v))/(k);
S = (k+mu)*(d+r+mu-b*w*v)*N/(beta*k);
I = b*w*N/(T) - (sigma+mu)*(k+mu)*(d+r+mu-b*w*v)*N/(beta*k*T);
L = ((d+r+mu-b*w*v)/(k))*I;
P = b*(1-w)*N/(mu) + sigma*(k+mu)*(d+r+mu-b*w*v)*N/(beta*mu*k) + (r/mu)*I;

display(num2str([S L I P]));

HBV = (L+I)/(S+L+I+P);

display(HBV);

Nt = X2(:,1)+X2(:,2)+X2(:,3)+X2(:,4)+X2(:,5);

%% Plotting

myworkingdir = cd; % access current working directory
plotsdir = 'plots_Ph'; % create 'plots' folder in the current directory
cd(plotsdir); % go to the plots folder directory

% counting the number of files with extension name '*.fig'
filecount=numel(dir('*.fig')); 

% create name for the figure with the counter and corresponding extension 
% names (e.g., plot#.fig, plot#.png)
mynewfigname1 = strcat('HepaB_Incidence_beta_k_timedep_2019',int2str(filecount+1),'.fig');
mynewpngname1 = strcat('HepaB_Incidence_beta_k_timedep_2019',int2str(filecount+1),'.png');


% create figure name (plot#)
plotname1 = strcat('HepaB_Incidence',int2str(filecount+1));
plotname2 = strcat('HepaB_Incidence',int2str(filecount+1));

xx1 = 2017:0.05:2019;
yy1 = spline(t_data(1:3),I1,xx1);
xx2 = 2020:0.05:2023;
yy2 = spline(t_data(4:end),I2,xx2);

figure(1); clf; hold on;
plot(xx1,yy1,'Color','b','LineWidth',2);
plot(t_data,I_data,'rs','LineWidth',3,'MarkerSize',8);
plot(xx2,yy2,'Color','b','LineWidth',2);
set(gca,'ytick',0:20000:90000,'fontsize',10)
legend('model','data');
xlabel('year','FontSize',16); ylabel('HBV Incidence','FontSize',16);
saveas(gcf,mynewfigname1)
saveas(gcf,mynewpngname1)


