function [Cbar]=Solve_Global_Damp(obj, node, C, U)
% SOLVE_GLOBAL_DAMP Assembles the global damping matrix Cbar for bar elements.

    nodalCoordinates = node.coordinates_mat;
    N_node = size(nodalCoordinates, 1);
    Cbar = zeros(3*N_node,3*N_node);
    
end
