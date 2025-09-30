function [Tbar,Kbar,Cbar]=Solve_FK(obj,node,U)
    
        [Ex]=Solve_Strain(obj,node,U);
        [Sx,C]=Solve_Stress(obj,Ex);
        [Tbar]=Solve_Global_Force(obj,node,U,Sx);
        [Kbar]=Solve_Global_Stiff(obj,node,U,Sx,C);
        [Cbar]=Solve_Global_Damp(obj, node, C, U);

        obj.strain_current_vec=Ex;
        obj.energy_current_vec=1/2*obj.E_vec.*obj.A_vec.*obj.L0_vec.*Ex.*Ex;

end