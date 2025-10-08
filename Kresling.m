%% Initialize the solver
% Khởi tạo bộ giải
clear all;
clc;
close all;

tic % Bắt đầu đếm thời gian thực hiện chương trình

%% Define the Geometry of origami
% Định nghĩa Hình học của mô hình origami
% Tại đây chúng ta tạo hình học của origami

% Định nghĩa tọa độ nút trước khi chia lưới
R=40*10^(-3); % Bán kính của hình trụ cơ bản (đơn vị: mét)
H=50*10^(-3); % Chiều cao của một tầng (đơn vị: mét)
theta=30/180*pi; % Góc xoắn (đơn vị: radian)
N=6; % Số mặt của đa giác đáy (6 cạnh)
m=1; % Số tầng của cấu trúc

% Tham số độ cứng của cấu trúc
sprStiff=0.000001; % Độ cứng của lò xo xoắn 0.000001

barE=1*10^9; % Young's modulus (Mô-đun Young)
KLR = 900; % Khối lượng riêng kg/m3 900
panelv=0.3; % Hệ số Poisson của các mặt phẳng
panelt=2*10^-3; % Độ dày của các mặt phẳng

% Tìm diện tích của các thanh (bars)
alpha=2*pi/N;
node1=[R*cos(theta),R*sin(theta),0];
node2=[R*cos(alpha+theta),R*sin(alpha+theta),0];
node3=[R*cos(alpha+theta*2),R*sin(alpha+theta*2),H];
v1=node2-node1;
v2=node3-node1;
TriangleArea=norm(cross(v2,v1))/2; %tìm diện tích mặt
TriangleLength=norm(v1)+norm(v2)+norm(node3-node2); %tìm chu vi mặt

barA=TriangleArea*panelt/TriangleLength*2/(1-panelv);
% Tổng diện tích của thanh được tính toán bằng cách bảo toàn tổng thể tích của mặt phẳng,
% giúp bảo toàn tổng năng lượng đàn hồi dưới biến dạng đồng nhất.
% Công thức này được lấy từ bài báo Merlin2 của Liu.


% Khởi tạo các phần tử
bar=Vec_Elements_Bars; % Đối tượng để lưu trữ các thanh (bars)
rotSpr=Vec_Elements_RotSprings_4N; % Đối tượng để lưu trữ các lò xo xoắn
node=Elements_Nodes; % Đối tượng để lưu trữ các nút (nodes)

%% Geometry of Kresling origami
% Hình học của origami Kresling

% Tạo tọa độ cho tất cả các nút
for i=1:m+1
    for j=1:N
        node.coordinates_mat=[node.coordinates_mat;
            R*cos(j*alpha+i*theta),R*sin(j*alpha+i*theta),(i-1)*H];
    end
end

% Thêm các thanh ngang
for i=1:m+1
    for j=1:N
        if j ~=N
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i-1)*N+j+1];
        else
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i-1)*N+1];
        end
    end
end

% Thêm các thanh chéo (diagonal bars)
for i=1:m
    for j=1:N
        if j ~=N
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+j];
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+j+1];
        else
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+j];
            bar.node_ij_mat=[bar.node_ij_mat;
                (i-1)*N+j,(i)*N+1];
        end
    end
end

% Khởi tạo các tham số khác
barNum=size(bar.node_ij_mat);
barNum=barNum(1);
bar.A_vec=barA*ones(barNum,1); % Gán diện tích tiết diện cho tất cả các thanh
bar.E_vec=barE*ones(barNum,1); % Gán mô-đun Young cho tất cả các thanh

%% Set up the rotational springs
% Thiết lập các lò xo xoắn
% Lò xo xoắn chéo (diagonal)
for i=1:m
    for j=1:N
        if j ==1
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+N,(i-1)*N+j,(i)*N+j,(i)*N+j+1];
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i)*N+j,(i-1)*N+j,(i)*N+j+1,(i-1)*N+j+1];
        elseif j~=N
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j-1,(i-1)*N+j,(i)*N+j,(i)*N+j+1];
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i)*N+j,(i-1)*N+j,(i)*N+j+1,(i-1)*N+j+1];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j-1,(i-1)*N+j,(i)*N+j,(i)*N+1];
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i)*N+j,(i-1)*N+j,(i)*N+1,(i-1)*N+1];
        end
    end
end

% Lò xo xoắn giữa các tầng
for i=1:m-1
    for j=1:N
        if j~=N
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j,i*N+j,i*N+j+1,(i+1)*N+j+1];
        else
            rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
                (i-1)*N+j,i*N+j,i*N+1,(i+1)*N+1];
        end
    end
end
rotSpr.rot_spr_K_vec=sprStiff*ones(m*(2*N)+N*(m-1),1); % Gán độ cứng cho tất cả các lò xo

%% Initialize assembly
% Khởi tạo lắp ráp
assembly=Assembly_Origami(); % Tạo đối tượng lắp ráp
assembly.node=node; % Gán các nút
assembly.bar=bar; % Gán các thanh
assembly.rotSpr=rotSpr; % Gán các lò xo xoắn

assembly.Initialize_Assembly() % Khởi tạo lắp ráp để chuẩn bị cho việc giải

%% Plot for investigation
% Vẽ đồ thị để kiểm tra
plots=Plot_Origami();
plots.displayRange=0.1;
plots.displayRangeRatio=1;
plots.assembly=assembly;

% plots.Plot_Shape_NodeNumber() % Vẽ hình dạng và hiển thị số nút
% plots.Plot_Shape_BarNumber() % Vẽ hình dạng và hiển thị số thanh
% plots.Plot_Shape_SprNumber() % Vẽ hình dạng và hiển thị số lò xo xoắn

% thiết lập các mặt phẳng để vẽ
panelNum=1;
for i=1:m
    for j=1:N
        if j ~=N
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i-1)*N+j+1,(i)*N+j+1];
            panelNum=panelNum+1;
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i)*N+j,(i)*N+j+1];
            panelNum=panelNum+1;
        else
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i-1)*N+1,(i)*N+1];
            panelNum=panelNum+1;
            plots.panelConnection{panelNum}=[
                (i-1)*N+j,(i)*N+j,(i)*N+1];
            panelNum=panelNum+1;
        end
    end
end
%%%%%%%%%%%%%%%
% Tạo ma trận khối lượng
node.mass_vec= (R * sin(pi/N) * sqrt( H^2 + 4*R^2 * (sin(theta/2))^2 * (sin(2*pi/N))^2 )) * panelt * KLR;
        M = node.FindMassMat;

%% Setup the loading controller
% Thiết lập bộ điều khiển tải
dc=Solver_DC;
dc.assembly=assembly;

dc.supp = [(1:N)'  ones(N,3)];


dc.dt = 0.00001;

force=120; % Lực tác dụng

dc.selectedRefDisp=[6*m+1,3]; % Chọn nút tham chiếu để đo biến dạng

dc.load=[6*m+1,0,0,-force; % Áp dụng lực nén lên các nút trên đỉnh
         6*m+2,0,0,-force;
         6*m+3,0,0,-force;
         6*m+4,0,0,-force;
         6*m+5,0,0,-force;
         6*m+6,0,0,-force];

dc.increStep=1500; % Số bước tăng tải
dc.tol=10^-5; % Sai số chấp nhận được
dc.iterMax=10000; % Số lần lặp tối đa

% [Uhis, Vhis, Ahis] = dc.Solve_byRK4(); % Chạy mô phỏng và lưu trữ lịch sử biến dạng
[Uhis, Vhis, Ahis] = dc.Solve_byNewmark_beta();
% [Uhis] = dc.Solve();
Uhis;

% plots.fileName='Kresling.gif';
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:))); % Vẽ hình dạng cuối cùng sau biến dạng
plots.Plot_DeformedHis(Uhis); % Vẽ lịch sử biến dạng (đã bị chú thích)

time = (0:dc.increStep-1) * dc.dt;  % vector thời gian
NodeNum = size(Uhis,2);

% Tính vận tốc tổng hợp
h = N+1;
Vmag = sqrt( squeeze(Vhis(:,h,1)).^2 + squeeze(Vhis(:,h,2)).^2 + squeeze(Vhis(:,h,3)).^2 );
% Vmag = squeeze(Vhis(:,h,3));

% Vẽ đồ thị
figure;
plot(time, Vmag, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Velocity Magnitude [m/s]');
title('Velocity Magnitude of Node');
grid on;



% % %% Find the reaction force and loading results
% % % Tìm lực phản ứng và kết quả tải
% % 
% % forceHis=zeros(dc.increStep,1);
% % UrefHis=zeros(dc.increStep,1);
% % 
% % for i=1:dc.increStep
% %     [F,K]=assembly.Solve_FK(squeeze(Uhis(i,:,:)));%F: véc-tơ lực nút nội bộ / lực mất cân bằng tương ứng với U;K: ma trận tiếp tuyến.
% %     UrefHis(i)=Uhis(i,dc.selectedRefDisp(1),dc.selectedRefDisp(2));%chuyển vị Z của nút tham chiếu ở bước i.
% %     forceHis(i)=F(dc.selectedRefDisp(2)+(dc.selectedRefDisp(1))*3);
% % end
% % 
% % figure
% % plot(-[0;UrefHis],[0;forceHis]) % Vẽ biểu đồ lực-biến dạng
% % xlabel('Z deformation of top node (m)') % Nhãn trục X
% % ylabel('Applied Force (N)') % Nhãn trục Y
% % 
% % toc % Dừng đếm thời gian
% 
% % Các dòng mã dưới đây được chú thích (comment) để tính toán và vẽ biểu đồ năng lượng biến dạng,
% % nhưng không được thực thi trong phiên bản hiện tại.
% % EnergyHis=zeros(dc.increStep,2);
% % 
% % for i=1:dc.increStep
% %     % rotational spring element
% %     rotSpr.SolveFK(node,squeeze(Uhis(i,:,:)));
% %     EnergyHis(i,1)=sum(rotSpr.currentStrainEnergy_Vec);
% % 
% %     % bar element
% %     bar.SolveFK(node,squeeze(Uhis(i,:,:)));
% %     EnergyHis(i,2)=sum(bar.currentStrainEnergy_Vec);
% % 
% % end
% % 
% % figure
% % hold on
% % plot(UrefHis,EnergyHis(:,1))
% % plot(UrefHis,EnergyHis(:,2))
% % xlabel('Z deformation of top node (m)') 
% % ylabel('Strain Energy (J)')
% % legend('rotational springs','bars')