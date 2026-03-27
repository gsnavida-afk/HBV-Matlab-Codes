function X = model_PRCC(t,X0,b,mu,d,beta,epsilon,k,r1,r2,r3,w,p,q,v)

S = X0(1);
L = X0(2);
I = X0(3);
C = X0(4);
R = X0(5);

N = S + L + I + C + R;

dS = b*w*(N - v*C) - beta*S./N*(I+epsilon*C) - (mu + r3)*S;
dL = beta*S./N*(I+epsilon*C) - (k+mu)*L;
dI = q*k*L - (mu+r1)*I;
dC = b*w*v*C + p*r1*I + (1-q)*k*L - (mu+d+r2)*C;
dR = b*(1-w)*N + (1-p)*r1*I + r2*C + r3*S - mu*R;

X = [dS;dL;dI;dC;dR];