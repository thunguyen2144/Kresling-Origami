%% Node Elements
% Here, we treat the node also as an element. This can be useful for
% dynamic simulation, where these nodes are like a lumped mass element. 
% The primary usage of node element is to store nodal cooridnation. 

classdef Elements_Nodes < handle

    properties
        % Nodal corridnates stored as a matrix (N*3)
        % Total number of node (N)
        coordinates_mat   % da co

        % Mass of nodes, stored as a vector (N*1)
        mass_vec  % co

        % Current deformation, stored as a Matrix (N*3)
        current_U_mat

        % Current external force (N*3)
        current_ext_force_mat
    end

    methods   
        % Find the mass matrix of the system
        M=FindMassMat(obj)
    end
end