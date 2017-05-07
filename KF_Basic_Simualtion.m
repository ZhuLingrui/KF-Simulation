N = 100;
w = randn(1,N);
X(1)  = 0;
A = 1;
for k = 2:N
    X(k) = A*X(k-1)+w(k-1);   % describe the proceess
end

V = randn(1,N);
q1 = std(V);
Rvv = q1.^2;
q2 = std(X);
Rxx = q2.^2;
q3 = std(w);
Rww = q3.^2;

H = 0.8;
Z = H*X+w;                 % decribe the measurement

X_hat_priori = zeros(1,N);   % Pre-allocation
P_priori = zeros(1,N);
K = zeros(1,N);
X_hat_post = zeros(1,N);
P_post = zeros(1,N);

for t = 2:N
    X_hat_priori(t) = A*X_hat_post(t-1);          % Time update
    P_priori(t) = A.^2*P_post(t-1) + Rww;
    K(t) = H*P_priori(t)/(H.^2*P_priori(t)+Rvv);    % measurement update
    X_hat_post(t) = X_hat_priori(t) + K(t)*(Z(t)-H*X_hat_priori(t));
    P_post(t) = (1-K(t)*H)*P_priori(t);
end

t = 1:N;
plot(t,X,'r',t,Z,'b',t,X_hat_post,'g');
title('KF Simulation');
legend('real state','measurement','estimation');

