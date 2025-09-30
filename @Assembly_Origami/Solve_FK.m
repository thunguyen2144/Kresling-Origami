function [T,K,C]=Solve_FK(obj,U)

    [Tbar,Kbar,Cbar]=obj.bar.Solve_FK(obj.node,U);
    [Trs,Krs,Cspr]=obj.rotSpr.Solve_FK(obj.node,U);
    
    T=Tbar+Trs;
    K = Kbar + Krs;
    C=Cbar+Cspr;
end