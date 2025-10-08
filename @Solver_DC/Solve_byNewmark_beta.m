%% Solve dynamic equation M*ddU + C*Udot + K*U = F(t) by Newmark_beta
function [Uhis, Vhis, Ahis] = Solve_byNewmark_beta(obj)

    % --- Parameters ---
    increStep   = obj.increStep;   % số bước thời gian
    dt          = obj.dt;          % bước thời gian
    assembly    = obj.assembly;
    supp        = obj.supp;
    load        = obj.load;
    
    % --- Node / DOF info ---
    NodeCount = size(assembly.node.coordinates_mat,1);
    ndof_full = NodeCount * 3;
    allDOFs = (1:ndof_full)';
    
    % --- Mass matrix (full) ---
    M = assembly.node.FindMassMat();

    % --- Assemble global external load vector (full) ---
    F_full = zeros(ndof_full,1);
    for i=1:size(load,1)
        n = load(i,1);
        idx = (n-1)*3 + (1:3);
        F_full(idx) = F_full(idx) + load(i,2:4)';
    end

    % --- Identify support (fixed) DOFs from obj.supp ---
    fixDofs = [];
    for i = 1:size(supp,1)
        nodeIdx = supp(i,1);
        if supp(i,2)==1, fixDofs(end+1) = 3*(nodeIdx-1)+1; end
        if supp(i,3)==1, fixDofs(end+1) = 3*(nodeIdx-1)+2; end
        if supp(i,4)==1, fixDofs(end+1) = 3*(nodeIdx-1)+3; end
    end
    


    freeDOFs = setdiff(allDOFs, fixDofs);
    % ndof_red = numel(freeDOFs);

    % --- Initial config (full) from node.current_U_mat ---
    U0_mat = assembly.node.current_U_mat; % NodeCount x 3
    U0_full = reshape(U0_mat.', [], 1);   % ndof_full x 1
    V0_full = zeros(ndof_full,1);

        [T,~, C] = assembly.Solve_FK(U0_mat);
    
    % --- Reduce M,C,K,F,T to free DOFs ---
    M_red = M(freeDOFs, freeDOFs);
    C_red = C(freeDOFs, freeDOFs);   
    % K_red = K_full(freeDOFs, freeDOFs);
    
    F_red = F_full(freeDOFs);
    T_red = T(freeDOFs);

    % --- Initial reduced fields ---
    U = U0_full(freeDOFs);
    V = V0_full(freeDOFs);
    A = M_red \ (F_red - T_red - C_red * V); % initial acceleration

    % --- Newmark parameters (constant average) ---
    beta = 1/4;
    gamma = 1/2;
    a0 = 1/(beta*dt^2);
    a1 = gamma/(beta*dt);
    a2 = 1/(beta*dt);
    a3 = 1/(2*beta) - 1;
    a4 = gamma/beta - 1;
    a5 = dt*(gamma/(2*beta) - 1);

    % --- History containers (full shape NodeCount x 3) ---
    Uhis = zeros(increStep, NodeCount, 3);
    Vhis = zeros(increStep, NodeCount, 3);
    Ahis = zeros(increStep, NodeCount, 3);

    % store initial step
    U_full_current = zeros(ndof_full,1);
    V_full_current = zeros(ndof_full,1);
    A_full_current = zeros(ndof_full,1);

    U_full_current(freeDOFs) = U;
    V_full_current(freeDOFs) = V;
    A_full_current(freeDOFs) = A;

    Uhis(1,:,:) = reshape(U_full_current,3,[]).';
    Vhis(1,:,:) = reshape(V_full_current,3,[]).';
    Ahis(1,:,:) = reshape(A_full_current,3,[]).';

    % --- Time loop (Modified Newton: update K,C each step but no inner Newton iterations) ---
    for step = 2:increStep

        % if external loads time-dependent, update F_full here (currently constant)
        Fext_red = F_red;

        % --- Re-evaluate internal matrices at current config (full) ---
        Umat_forFK = reshape(U_full_current,3,[]).'; % NodeCount x 3
        
            [T, K_full, C_full] = assembly.Solve_FK(Umat_forFK);

        % reduced versions
        T_red = T(freeDOFs);
        K_red = K_full(freeDOFs, freeDOFs);
        C_red = C_full(freeDOFs, freeDOFs);

        % --- Effective stiffness & RHS (Newmark implicit) ---
        K_eff = K_red + a1 * C_red + a0 * M_red;
        RHS = Fext_red  + M_red*(a0*U + a2*V + a3*A) + C_red*(a1*U + a4*V + a5*A);

        % Solve for U_new (reduced)
        % Use backslash (stable)
        % denta_U = K_eff \ RHS;
        % U_new = denta_U + U;
        U_new = K_eff \ RHS;

        % update acceleration and velocity (reduced)
        A_new = a0*(U_new - U) - a2*V - a3*A;
        V_new = V + dt*((1 - gamma)*A + gamma*A_new);

        % map to full
        U_full_current = zeros(ndof_full,1);
        V_full_current = zeros(ndof_full,1);
        A_full_current = zeros(ndof_full,1);
        
        U_full_current(freeDOFs) = U_new;
        V_full_current(freeDOFs) = V_new;
        A_full_current(freeDOFs) = A_new;

        % store history
        Uhis(step,:,:) = reshape(U_full_current,3,[]).';
        Vhis(step,:,:) = reshape(V_full_current,3,[]).';
        Ahis(step,:,:) = reshape(A_full_current,3,[]).';

        % update reduced variables for next step 
         U = U_new; V = V_new; A = A_new;

    end

end
