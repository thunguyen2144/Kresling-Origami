function  [Sx,Cx]= Solve_Stress(obj,Ex)

    Sx=obj.E_vec.*Ex;
    Cx=obj.E_vec;

end