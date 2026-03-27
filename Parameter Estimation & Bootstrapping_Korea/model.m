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