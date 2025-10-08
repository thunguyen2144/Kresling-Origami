clc; clear; close all;

%% Tham số hệ 1 bậc tự do
M = 10; 
C = 1000; 
K = 1e6; 

dt = 0.0001; 
t_total = 0.1; 
t = 0:dt:t_total;
increStep = length(t);

% Điều kiện ban đầu
U0 = 0; 
V0 = 0;

% Lực ngoài (hằng)
F = 1;  % hoặc thử: F = sin(t(i));

%% Khởi tạo
U = U0; 
V = V0;

Uhis = zeros(increStep,1);
Vhis = zeros(increStep,1);
Ahis = zeros(increStep,1);

Uhis(1) = U;
Vhis(1) = V;
Ahis(1) = (F - C*V - K*U)/M;

%% --- Hàm đạo hàm hệ ---
% y = [U; V]
f = @(t, y) [ y(2);
              (F - C*y(2) - K*y(1)) / M ];

%% --- Vòng lặp RK4 ---
for i = 2:increStep
    y = [U; V];

    k1 = f(t(i), y);
    k2 = f(t(i) + dt/2, y + dt/2 * k1);
    k3 = f(t(i) + dt/2, y + dt/2 * k2);
    k4 = f(t(i) + dt,   y + dt * k3);

    y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4);

    % Cập nhật
    U = y_new(1);
    V = y_new(2);
    A = (F - C*V - K*U) / M;

    % Lưu kết quả
    Uhis(i) = U;
    Vhis(i) = V;
    Ahis(i) = A;
end

%% --- Vẽ kết quả ---
figure;
plot(t, Uhis, 'b-', 'LineWidth', 1.5);
title('Chuyển vị u(t) theo RK4');
xlabel('Thời gian (s)'); ylabel('u(t)');
grid on;

figure;
plot(t, Vhis, 'r-', 'LineWidth', 1.5);
title('Vận tốc v(t) theo RK4');
xlabel('Thời gian (s)'); ylabel('v(t)');
grid on;

figure;
plot(t, Ahis, 'k-', 'LineWidth', 1.5);
title('Gia tốc a(t) theo RK4');
xlabel('Thời gian (s)'); ylabel('a(t)');
grid on;
