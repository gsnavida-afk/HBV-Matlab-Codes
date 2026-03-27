function X = model_PRCC(t,X0,b,mu,d,beta,k,r,sigma,w,v,dummy)

S = X0(1);
L = X0(2);
I = X0(3);
R = X0(4);

N = S + L + I + R;


dS = b*w*(N - v*I) - beta*S./N*(I) - (mu + sigma)*S;
dL = beta*S./N*(I) - (k+mu)*L;
dI = b*w*v*I + k*L - (mu+d+r)*I; 
dR = b*(1-w)*N + r*I + sigma*S - mu*R;

X = [dS; dL; dI; dR];
end