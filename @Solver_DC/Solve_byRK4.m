%% Solve dynamic equation M*ddU + C*Udot + K*U = F(t) by RK4
function [Uhis, Vhis, Ahis] = Solve_byRK4(obj)

    % Parameters
    increStep   = obj.increStep;   % số bước thời gian
    dt          = obj.dt;          % bước thời gian
    assembly    = obj.assembly;
    supp        = obj.supp;
    load        = obj.load;

    % Initial condition (NodeNum x 3)
    U0 = assembly.node.current_U_mat;
    NodeNum = size(U0,1);

    % Mass matrix (full)
    M = assembly.node.FindMassMat();   % Kích thước ndof x ndof
    ndof_full = 3 * NodeNum;
    allDOFs = (1:ndof_full)';

    % --- Assemble load vector (full) ---
    F_full = zeros(ndof_full,1);
    for i=1:size(load,1)
        n = load(i,1);
        idx = (n-1)*3 + (1:3);
        F_full(idx) = F_full(idx) + load(i,2:4)';
    end

    % === Xử lý điều kiện biên (supp): xác định fixDofs và freeDOFs ===
    fixDofs = [];
    for i = 1:size(supp,1)
        nodeIdx = supp(i,1);
        if supp(i,2)==1, fixDofs(end+1) = 3*(nodeIdx-1)+1; end
        if supp(i,3)==1, fixDofs(end+1) = 3*(nodeIdx-1)+2; end
        if supp(i,4)==1, fixDofs(end+1) = 3*(nodeIdx-1)+3; end
    end
    fixDofs = unique(fixDofs);
    freeDOFs = setdiff(allDOFs, fixDofs);
    ndof_red = numel(freeDOFs);

    % --- Reduce mass and force for solver ---
    M_red = M(freeDOFs, freeDOFs);
    F_red = F_full(freeDOFs);

    % Initial reduced vectors
    U_full_init = reshape(U0.',[],1);    % ndof_full x 1 (column vector)
    U_red = U_full_init(freeDOFs);       % reduced initial displacement
    V_full_init = zeros(ndof_full,1);
    V_red = V_full_init(freeDOFs);

    % Storage (full shape)
    Uhis = zeros(increStep, NodeNum, 3);
    Vhis = zeros(increStep, NodeNum, 3);
    Ahis = zeros(increStep, NodeNum, 3);

    % Save initial state (step 1)
    U_full_current = zeros(ndof_full,1);
    U_full_current(freeDOFs) = U_red;
    V_full_current = zeros(ndof_full,1);
    V_full_current(freeDOFs) = V_red;

    % compute initial acceleration (need T,C at U0)
    [Tvec_init,~,C_init] = assembly.Solve_FK(reshape(U_full_current,3,[]).'); % T: internal nodal force vector (N?), K,C mats
    
    % Remove fixed DOFs effect (we treat supports as zero disp)
    T_red_init = Tvec_init(freeDOFs);
    C_red_init = C_init(freeDOFs, freeDOFs);
    Avec_red = M_red \ (F_red - T_red_init - C_red_init * V_red);

    % Store initial
    Uhis(1,:,:) = reshape(U_full_current,3,[]).';
    Vhis(1,:,:) = reshape(V_full_current,3,[]).';
    Ahis(1,:,:) = reshape(zeros(ndof_full,1),3,[]).'; % we'll fill from reduced
    A_full_current = zeros(ndof_full,1);
    A_full_current(freeDOFs) = Avec_red;
    Ahis(1,:,:) = reshape(A_full_current,3,[]).';

    % --- Define RHS function for reduced system Y_red = [U_red; V_red] ---
    function dY = f_reduced(t, Yred)
        % Unpack reduced U,V
        Ured = Yred(1:ndof_red);
        Vred = Yred(ndof_red+1:end);

        % Reconstruct full U for assembly.Solve_FK
        Ufull_temp = zeros(ndof_full,1);
        Ufull_temp(freeDOFs) = Ured;
        % Reshape to NodeNum x 3 for Solve_FK
        Umat = reshape(Ufull_temp, 3, []).';
        [Tmat,Kfull,Cfull] = assembly.Solve_FK(Umat);

        % Reduce T,K,C to free DOFs
        Tred = Tmat(freeDOFs);
        % Kred = Kfull(freeDOFs, freeDOFs);
        Cred = Cfull(freeDOFs, freeDOFs);

        % Compute acceleration (reduced)
        rhs_red = F_red - Tred - Cred * Vred;
        Ared = M_red \ rhs_red;

        % dUdt = V, dVdt = A
        dY = [Vred; Ared];
    end

    % --- Time integration RK4 on reduced system ---
    t = 0;
    Yred = [U_red; V_red];

    for step = 2:increStep
        % RK4 stages
        k1 = f_reduced(t, Yred);
        k2 = f_reduced(t + dt/2, Yred + dt/2 * k1);
        k3 = f_reduced(t + dt/2, Yred + dt/2 * k2);
        k4 = f_reduced(t + dt,   Yred + dt * k3);

        Yred_next = Yred + dt/6 * (k1 + 2*k2 + 2*k3 + k4);

        % Update
        Yred = Yred_next;
        t = t + dt;

        % Extract reduced U,V,A (compute A by one f_reduced call to ensure consistency)
        Ured = Yred(1:ndof_red);
        Vred = Yred(ndof_red+1:end);
        dY = f_reduced(t, Yred);
        Ared = dY(ndof_red+1:end);

        % Map back to full for storage
        % U_full_current = zeros(ndof_full,1);
        % V_full_current = zeros(ndof_full,1);
        % A_full_current = zeros(ndof_full,1);

        U_full_current(freeDOFs) = Ured;
        V_full_current(freeDOFs) = Vred;
        A_full_current(freeDOFs) = Ared;

        Uhis(step,:,:) = reshape(U_full_current,3,[]).';
        Vhis(step,:,:) = reshape(V_full_current,3,[]).';
        Ahis(step,:,:) = reshape(A_full_current,3,[]).';
    end


end
