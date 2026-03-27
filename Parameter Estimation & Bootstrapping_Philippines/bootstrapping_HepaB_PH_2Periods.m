clc;clear; close all
load IncidenceHepaBPh1.mat %Hepa B prevalence data from 2010-2023

t_data = IncidenceHepaBPh1(:,1);
I_data = IncidenceHepaBPh1(:,2);

%% Parameter Estimation on original data for the First Period (2017-2019)

% We estimate the values of the parameters 
% beta, r

ini_param=[2.6797 0.025]; LB=[0 0]; UB=[5 1];
options = optimoptions("lsqcurvefit",'StepTolerance',1e-15,'MaxFunctionEvaluations',1e+5);

% Original data
xvalues_orig1 = lsqcurvefit(@model_fit1,ini_param,t_data(1:3),I_data(1:3),LB,UB,options);
beta_orig1 =xvalues_orig1(1);
r_orig1 =xvalues_orig1(2);

display(xvalues_orig1);

%% Best Fit Model for Original Data for the First Period (2017-2019)

% Other Parameter Values
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

X0 = [S0; L0; I0; R0; I02];
options = [];

[~,Xsamp2] = ode45(@modelaux,t_data(1:3),X0,[],b,mu,d,beta_orig1,k,r_orig1,sigma,w,v);
I2 = [Xsamp2(1,5); diff(Xsamp2(:,5))];

%% Generate Synthetic Data for the First Period (2017-2019) from the best-fit model
n = 10000; %number of synthetic datasets to be generated

beta1 = zeros([1 n]);
r1 = zeros([1 n]);

D_synthetic = zeros([3,n]);
for i = 1:n
    %D_synthetic(:,i) = normrnd(I_data,1); %generate synthetic data using normal random numbers
    %D_synthetic(:,i) = nbinrnd(I_data,0.5); %generate synthetic data using negative binomial random numbers
    D_synthetic(:,i) = poissrnd(I2); %generate synthetic data using Poisson random numbers
end

% Parameter Estimation for the Synthetic Datasets
for i = 1:n
    xvalues1 = lsqcurvefit(@model_fit1,ini_param,t_data(1:3),D_synthetic(:,i),LB,UB,options);
    beta1(i) = xvalues1(1);
    r1(i) = xvalues1(2);
end


figure(3); clf; hold on; 
for i=1:n

    [~,X1samp] = ode45(@modelaux,t_data(1:3),X0,[],b,mu,d,beta1(i),k,r1(i),sigma,w,v);
    I1 = [X1samp(1,5); diff(X1samp(:,5))];
    
    plot(t_data(1:3),I1,'c-')
end


%% Parameter Estimation on original data and each synthetic data set for the Second Period (2020-2023)

ini_param=[2.6797 0.025]; LB=[0 0]; UB=[5 1];
options = optimoptions("lsqcurvefit",'StepTolerance',1e-15,'MaxFunctionEvaluations',1e+5);

% Original data
xvalues_orig2 = lsqcurvefit(@model_fit2,ini_param,t_data(4:end),I_data(4:end),LB,UB,options);
beta_orig2 =xvalues_orig2(1);
r_orig2 =xvalues_orig2(2);

display(xvalues_orig1);

%% Best Fit Model for Original Data for the Second Period (2020-2023)

N0 = 106700000;
I0 = 40754;
L0 = I0/4;
R0 = 0.3*N0;
I02 = 40754;
S0 = N0-L0-I0-R0;

X0 = [S0; L0; I0; R0; I02];
options = [];

tspan=(2020:1:2030);

[~,Xsamp2] = ode45(@modelaux,tspan,X0,[],b,mu,d,beta_orig2,k,r_orig2,sigma,w,v);
I4 = [Xsamp2(1,5); diff(Xsamp2(:,5))];


%% Generate Synthetic Data for the Second Period (2020-2023) from the best fit model

beta2 = zeros([1 n]);
r2 = zeros([1 n]);

E_synthetic = zeros([4,n]);
for i = 1:n
    %D_synthetic(:,i) = normrnd(I_data,1); %generate synthetic data using normal random numbers
    %D_synthetic(:,i) = nbinrnd(I_data,0.5); %generate synthetic data using negative binomial random numbers
    E_synthetic(:,i) = poissrnd(I4(1:4)); %generate synthetic data using Poisson random numbers
end

% Parameter estimation for Synthetic Data
for i = 1:n
    xvalues2 = lsqcurvefit(@model_fit2,ini_param,t_data(4:end),E_synthetic(:,i),LB,UB,options);
    beta2(i) = xvalues2(1);
    r2(i) = xvalues2(2);
end

figure(3); hold on; 
for i=1:n

    [~,X1samp] = ode45(@modelaux,tspan,X0,[],b,mu,d,beta2(i),k,r2(i),sigma,w,v);
    I3 = [X1samp(1,5); diff(X1samp(:,5))];
    
    plot(tspan,I3,'c-')
end

figure(3); hold on; 
plot(t_data(1:3),I2,'k')
plot(t_data,I_data,'rs','LineWidth',3,'MarkerSize',8)
plot(tspan,I4,'k')
xline(2023,'-.','Projection','DisplayName','Projection')
set(gca,'ytick',0:20000:110000,'fontsize',10)
xlim([2017 2030])
xlabel('Year') 
ylabel('HBV Incidence')


% mean, std and  95% CI from distribution of parameter estimates
param_beta1=[mean(beta1) std(beta1)  plims(beta1,0.025) plims(beta1,0.975)];
param_beta2=[mean(beta2) std(beta2)  plims(beta2,0.025) plims(beta2,0.975)];
param_r1=[mean(r1) std(r1) plims(r1,0.025) plims(r1,0.975)];
param_r2=[mean(r2) std(r2) plims(r2,0.025) plims(r2,0.975)];
%% Histogram of Range of Parameters
numBins = 20;
ylimadd = 100;
figure(1)
tiledlayout(2,2)

ax1 = nexttile;
hist1 = histogram(beta1,numBins,'EdgeColor','k',FaceColor='#80B3FF');
grid off;
xlabel('\beta_{A}', 'FontSize', 15);
mu = mean(beta1);
sigma = std(beta1);
xline(beta_orig1, 'Color', 'b', 'LineWidth', 2);
xline(mu, 'Color', 'g', 'LineWidth', 2);
xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
ylim([0, max(hist1.Values)+ylimadd]); % Give some headroom above the bars.
sMean1= sprintf('Mean = %.3e.  SD = %.3e', mu, sigma);
title(sMean1, 'FontSize', 10);

ax2 = nexttile;
hist2 = histogram(beta2,numBins,'EdgeColor','k',FaceColor='#80B3FF');
grid off;
xlabel('\beta_{B}', 'FontSize', 15);
mu = mean(beta2);
sigma = std(beta2);
xline(beta_orig2, 'Color', 'b', 'LineWidth', 2);
xline(mu, 'Color', 'g', 'LineWidth', 2);
xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
ylim([0, max(hist2.Values)+ylimadd]); % Give some headroom above the bars.
sMean1= sprintf('Mean = %.3e.  SD = %.3e', mu, sigma);
title(sMean1, 'FontSize', 10);


ax3 = nexttile;
hist3 = histogram(r1,numBins,'EdgeColor','k',FaceColor='#80B3FF');
grid off;
xlabel('\rho_{A}', 'FontSize', 15);
mu = mean(r1);
sigma = std(r1);
xline(r_orig1, 'Color', 'b', 'LineWidth', 2);
xline(mu, 'Color', 'g', 'LineWidth', 2);
xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
ylim([0, max(hist3.Values)+ylimadd]); % Give some headroom above the bars.
sMean1= sprintf('Mean = %.3e.  SD = %.3e', mu, sigma);
title(sMean1, 'FontSize', 10);

ax4 = nexttile;
hist4 = histogram(r2,numBins,'EdgeColor','k',FaceColor='#80B3FF');
grid off;
xlabel('\rho_{B}', 'FontSize', 15);
mu = mean(r2);
sigma = std(r2);
xline(r_orig2, 'Color', 'b', 'LineWidth', 2);
xline(mu, 'Color', 'g', 'LineWidth', 2);
xline(mu - sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
xline(mu + sigma, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
ylim([0, max(hist4.Values)+ylimadd]); % Give some headroom above the bars.
sMean1= sprintf('Mean = %.3e.  SD = %.3e', mu, sigma);
title(sMean1, 'FontSize', 10);






