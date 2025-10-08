clc; clear; close all;

%% Tham số hệ 1 bậc tự do
M = 10; 
C = 1000; 
K = 1000000; 

dt = 0.0001; 
t_total = 0.1; 
t = 0:dt:t_total;
increStep = length(t);

% Điều kiện ban đầu
U0 = 0; 
V0 = 0;

% Lực ngoài
% F = sin(t)';
F = 1;
% Tham số Newmark
beta  = 0.25;
gamma = 0.5;

%% Khởi tạo
U = U0; 
V = V0;
A = (F(1) - C*V - K*U)/M;   % gia tốc ban đầu

Uhis = zeros(increStep,1);
Vhis = zeros(increStep,1);
Ahis = zeros(increStep,1);

Uhis(1) = U;
Vhis(1) = V;
Ahis(1) = A;

%% Các hệ số Newmark
a0 = 1/(beta*dt^2);
a1 = gamma/(beta*dt);
a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;
a5 = dt*(gamma/(2*beta) - 1);

Keff = K + a0*M + a1*C
%% Time stepping
for i = 2:increStep
    % RHS
    RHS = F  + M*(a0*U + a2*V + a3*A) + C*(a1*U + a4*V + a5*A);
    % RHS = F + M*(a0*U + a2*V + a3*A) + C*(a1*U + a4*V + a5*A);
    % Giải U_{i+1}
    Unew = RHS / Keff;

    % Tính A_{i+1}, V_{i+1}
    Anew = a0*(Unew - U) - a2*V - a3*A;
    Vnew = V + dt*((1-gamma)*A + gamma*Anew);

    % Cập nhật
    U = Unew;
    V = Vnew;
    A = Anew;

    % Lưu kết quả
    Uhis(i) = U;
    Vhis(i) = V;
    Ahis(i) = A;
end

%% Vẽ kết quả
figure;
plot(t,Uhis,'b-','LineWidth',1.5); hold on;
title('Chuyển vị u(t) với F=1');
xlabel('Thời gian (s)'); ylabel('u(t)');
grid on;

figure;
plot(t,Vhis,'r-','LineWidth',1.5); hold on;
title('Vận tốc v(t)');
xlabel('Thời gian (s)'); ylabel('v(t)');
grid on;

figure;
plot(t,Ahis,'k-','LineWidth',1.5); hold on;
title('Gia tốc a(t)');
xlabel('Thời gian (s)'); ylabel('a(t)');
grid on;
