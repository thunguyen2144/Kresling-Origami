function [Cspr] = Solve_Global_Damp(obj, node, U)
    nodalCoordinates = node.coordinates_mat;
    N_node = size(nodalCoordinates, 1);
    Cspr = zeros(3*N_node,3*N_node);
    
end
