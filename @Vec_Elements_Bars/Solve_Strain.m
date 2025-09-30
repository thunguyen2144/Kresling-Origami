function [Ex]=Solve_Strain(obj,node,U)

    nodalCoordinates=node.coordinates_mat;
    barConnect=obj.node_ij_mat;
    barLength=obj.L0_vec;

    NodeIndex1=barConnect(:,1);
    NodeIndex2=barConnect(:,2);
    
    node1=nodalCoordinates(NodeIndex1,:);
    node2=nodalCoordinates(NodeIndex2,:);
    
    B1n=(1./(barLength.*barLength)).*[-(node2-node1) (node2-node1)];  
    
    iden=eye(3);
    idenMat=[iden -iden; -iden iden];
    
    Utemp=[U(NodeIndex1,:)';U(NodeIndex2,:)'];     
    B2Utemp=(1./(barLength.*barLength)).*(idenMat*Utemp)';  
    Ex=dot(B1n,Utemp',2)+0.5*dot(B2Utemp,Utemp',2);
    
end