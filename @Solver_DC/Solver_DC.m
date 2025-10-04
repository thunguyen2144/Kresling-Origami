%% Displacement Controlled Implicity Solver

classdef Solver_DC  < handle

    properties

        % assembly of the structure
        assembly

        dt

        % storing the support information
        supp
        
        % the applied load
        load
        % a sampel code is the following
        % loadForce=3;
        % load=[29,0,0,loadForce;
        %       30,0,0,loadForce;
        %       31,0,0,loadForce;
        %       32,0,0,loadForce;];
        
        % the total number of incremental steps
        increStep
        
        % the tolerance for each iteration
        tol=1*10^-5
        
        % the lambdaBar used to control the loading
        lambdaBar=1
        
        % the maximum allowed iteration number
        iterMax=30        
       
        % The selected reference displacement 
        selectedRefDisp        
        
        % the history of displacement field
        Uhis
                   
    end

    methods
        % Solve for the equilibrium results
        [Uhis, Vhis, Ahis] = Solve(obj);

    end



end