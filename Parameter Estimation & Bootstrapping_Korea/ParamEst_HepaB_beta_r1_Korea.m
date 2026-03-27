% Parameter Estimation

clc;clear; close all
load PrevalenceDataHepaB2017to2023.mat %Hepa B prevalence data from 2010-2023

t_data = PrevalenceDataHepaB2017to2023(:,1);
I_data = PrevalenceDataHepaB2017to2023(:,2);
C_data = PrevalenceDataHepaB2017to2023(:,3);
%% Optimization
init_par = [2.6797 2]; % initial value of birth rate
LB = [0 0];
UB = [5 5];
options = optimoptions("lsqcurvefit",'StepTolerance',1e-15,'MaxFunctionEvaluations',1e+6);
xopt = lsqcurvefit(@model_fit,init_par,t_data,[I_data C_data],LB,UB,options);

beta = xopt(1);
r1 = xopt(2);

display(xopt);
%% Best Fit Model

b = [2.7207e-02,2.269e-02]; % estimated
mu=[1/43.26650546, 1/41.6608984]; % estimated
d=0.001; 
epsilon=0.16;
k=4;
r2=0.025;
r3 = 0.01;
p=0.05;
q=0.1;
v=0.03;
w=0.05;

N0 = 51230704;
I0 = 6188;
L0 = I0/4;
C0 = 424359;
R0 = 0.3*N0;
S0 = N0-L0-I0-C0-R0;

init = [S0; L0; I0; C0; R0];

options = [];
tspan = linspace(2017,2023,1000);
[~,X] = ode45(@model,tspan,init,options,b,mu,d,beta,epsilon,k,r1,r2,r3,w,p,q,v);

Nt = X(:,1)+X(:,2)+X(:,3)+X(:,4)+X(:,5);

R01 = beta*b(1)*k*w/((k+mu(1))*(r1+mu(1))*(r3+mu(1)));
R02 = q + epsilon*(p*q*r1+(r1+mu(1))*(1-q))/(d+r2+mu(1)-b(1)*w*v);
R01st = R01*R02;

display(num2str(R01st));

R01 = beta*b(2)*k*w/((k+mu(2))*(r1+mu(2))*(r3+mu(2)));
R02 = q + epsilon*(p*q*r1+(r1+mu(2))*(1-q))/(d+r2+mu(2)-b(2)*w*v);
R02nd = R01*R02;

display(num2str(R02nd));


%% Plotting

myworkingdir = cd; % access current working directory
plotsdir = 'plots_Korea_2017to2023_edited'; % create 'plots' folder in the current directory
cd(plotsdir); % go to the plots folder directory

% counting the number of files with extension name '*.fig'
filecount=numel(dir('*.fig')); 

% create name for the figure with the counter and corresponding extension 
% names (e.g., plot#.fig, plot#.png)
mynewfigname1 = strcat('Acute_HepaB_beta_r1_5',int2str(filecount+1),'.fig');
mynewpngname1 = strcat('Acute_HepaB_beta_r1_5',int2str(filecount+1),'.png');
mynewfigname2 = strcat('Chronic_HepaB_beta_r1_5',int2str(filecount+1),'.fig');
mynewpngname2 = strcat('Chronic_HepaB_beta_r1_5',int2str(filecount+1),'.png');

% create figure name (plot#)
plotname1 = strcat('Acute_HepaB_beta_r1',int2str(filecount+1));
plotname2 = strcat('Chronic_HepaB_beta_r1',int2str(filecount+1));

figure(1); clf; hold on;
plot(t_data,I_data,'rs','LineWidth',3,'MarkerSize',8);
plot(tspan,X(:,3),'Color','b','LineWidth',2);
set(gca,'ytick',0:1000:6500,'fontsize',10)
xlabel('year','FontSize',16); ylabel('Acutely Infected','FontSize',16);
saveas(gcf,mynewfigname1)
saveas(gcf,mynewpngname1)

figure(2); clf; hold on;
plot(t_data,C_data,'rs','LineWidth',3,'MarkerSize',8);
plot(tspan,X(:,4),'Color','b','LineWidth',2);
set(gca,'ytick',0:10000:455000,'fontsize',10)
xlabel('year','FontSize',16); ylabel('Chronically Infected','FontSize',16);
saveas(gcf,mynewfigname2)
saveas(gcf,mynewpngname2)


%% Objective cost function
function Y = model_fit(x,year)
beta = x(1);
r1 = x(2);


b = [2.7207e-02,2.269e-02]; % estimated
mu=[1/43.26650546, 1/41.6608984]; % estimated
d=0.001; 
epsilon=0.16;
k=4;
r2=0.025;
r3 = 0.01;
p=0.05;
q=0.1;
v=0.03;
w=0.05;


N0 = 51230704;
I0 = 6188;
L0 = I0/4;
C0 = 424359;
R0 = 0.3*N0;
S0 = N0-L0-I0-C0-R0;

init = [S0; L0; I0; C0; R0];

options = [];
[~,HB_model] = ode45(@(t,x) ...
    model(t,x,b,mu,d,beta,epsilon,k,r1,r2,r3,w,p,q,v),year,init,options);
Y1 = HB_model(:,3);
Y2 = HB_model(:,4);
Y = [Y1 Y2];
end
%% Model
function X = model(t,X0,b,mu,d,beta,epsilon,k,r1,r2,r3,w,p,q,v)

S = X0(1);
L = X0(2);
I = X0(3);
C = X0(4);
R = X0(5);

N = S + L + I + C + R;

by = b(1);
if t > 2020
    by = b(2);
end

mu0 = mu(1);
if t > 2020
    mu0 = mu(2);
end

dS = by*w*(N - v*C) - beta*S./N*(I+epsilon*C) - (mu0 + r3)*S;
dL = beta*S./N*(I+epsilon*C) - (k+mu0)*L;
dI = q*k*L - (mu0+r1)*I; 
dC = by*w*v*C + p*r1*I + (1-q)*k*L - (mu0+d+r2)*C;
dR = by*(1-w)*N + (1-p)*r1*I + r2*C + r3*S - mu0*R;

X = [dS; dL; dI; dC; dR];
end