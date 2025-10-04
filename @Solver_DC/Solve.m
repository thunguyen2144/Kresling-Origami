%% Solve dynamic equation M*ddU + C*Udot + K*U = F(t) by RK4
function [Uhis, Vhis, Ahis] = Solve(obj)

    % Parameters
    increStep   = obj.increStep;   % số bước thời gian
    dt          = obj.dt;         % bước thời gian
    assembly    = obj.assembly;
    supp        = obj.supp;
    load        = obj.load;

    % Initial condition
    U = assembly.node.current_U_mat;   % (NodeNum x 3)
    V = zeros(size(U));                % vận tốc ban đầu = 0
    NodeNum = size(U,1);

    % Mass matrix
    M = assembly.node.FindMassMat();

    % Assemble load vector
    loadVec = zeros(3*NodeNum,1);
    for i=1:size(load,1)
        TempNodeNum = load(i,1);
        loadVec(TempNodeNum*3-2:TempNodeNum*3) = load(i,2:4);
    end

    % === Xử lý điều kiện biên (supp) ===
    fixDofs = [];
    for i = 1:size(supp,1)
        node = supp(i,1);
        if supp(i,2)==1, fixDofs(end+1) = 3*(node-1)+1; end
        if supp(i,3)==1, fixDofs(end+1) = 3*(node-1)+2; end
        if supp(i,4)==1, fixDofs(end+1) = 3*(node-1)+3; end
    end

    % Force BC
    loadVec(fixDofs) = 0;

    % Storage
    Uhis = zeros(increStep, NodeNum, 3);
    Vhis = zeros(increStep, NodeNum, 3);
    Ahis = zeros(increStep, NodeNum, 3);

    % Flatten U, V thành vector cột
    Uvec = reshape(U.',[],1);
    Vvec = reshape(V.',[],1);

    % Hàm f(t,Y) với Y = [U; V]
    function dY = f(t,Y)
        Uloc = Y(1:3*NodeNum);
        Vloc = Y(3*NodeNum+1:end);
        % Nội lực + độ cứng + ma trận cản
        [T,K,C] = assembly.Solve_FK(reshape(Uloc,3,[]).'); 
        [~,T] = Mod_K_For_Supp(K,supp,T);
        % Unbalance force
        Fext = loadVec;
        rhs = Fext - T - C*Vloc;
        A = M \ rhs; % gia tốc
        dY = [Vloc; A];
    end

    % Time integration RK4
    t = 0;
    for i=1:increStep
        Y = [Uvec; Vvec];

        k1 = f(t, Y);
        k2 = f(t + dt/2, Y + dt/2*k1);
        k3 = f(t + dt/2, Y + dt/2*k2);
        k4 = f(t + dt,   Y + dt*k3);

        Ynext = Y + dt/6*(k1 + 2*k2 + 2*k3 + k4);

        Uvec = Ynext(1:3*NodeNum);
        Vvec = Ynext(3*NodeNum+1:end);

        % Tính lại gia tốc tại bước này
        [T,K,C] = assembly.Solve_FK(reshape(Uvec,3,[]).');
        [~,T] = Mod_K_For_Supp(K,supp,T);
        rhs = loadVec - T - C*Vvec;
        Avec = M \ rhs;

        % === Gán lại điều kiện biên ===
        % Uvec(fixDofs) = 0;
        % Vvec(fixDofs) = 0;
        % Avec(fixDofs) = 0;

        % Lưu lại
        Uhis(i,:,:) = reshape(Uvec,3,[]).';
        Vhis(i,:,:) = reshape(Vvec,3,[]).';
        Ahis(i,:,:) = reshape(Avec,3,[]).';

        t = t + dt;
    end
end
