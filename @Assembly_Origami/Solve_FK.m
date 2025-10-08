function [T,K,C]=Solve_FK(obj,U)

    [Tbar,Kbar,Cbar]=obj.bar.Solve_FK(obj.node,U);
    [Trs,Krs,Cspr]=obj.rotSpr.Solve_FK(obj.node,U);
    
    T=Tbar+Trs;
    K = Kbar + Krs;
    C=Cbar+Cspr;

    % K = 50*eye(36,36);
    % C = 0.5*eye(36,36);
    % T = K*reshape(U.', [], 1)

end