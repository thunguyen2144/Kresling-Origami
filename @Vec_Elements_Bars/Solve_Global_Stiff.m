function [Kbar]=Solve_Global_Stiff(obj,node,U,Sx,C)
    
    nodalCoordinates=node.coordinates_mat;
    barConnect=obj.node_ij_mat;
    barLength=obj.L0_vec;
    barArea=obj.A_vec;

    A=size(nodalCoordinates);
    N_node=A(1); % Number of nodes
    
    A=size(C);
    N_bar=A(1); % Number of bars
    
    Kbar=zeros(3*N_node,3*N_node);
   
    NodeIndex1_temp=barConnect(:,1);
    NodeIndex2_temp=barConnect(:,2);
    node1_temp=nodalCoordinates(NodeIndex1_temp,:);
    node2_temp=nodalCoordinates(NodeIndex2_temp,:);
    
    barLength_square=barLength.*barLength;
    
    B1_temp=(1./barLength_square).*([-(node2_temp-node1_temp), (node2_temp-node1_temp)]);
    iden=eye(3);
    B2_pattern=[iden -iden; -iden iden];
    
    B2_temp=[1./barLength_square.*B2_pattern(1,:),...
        1./barLength_square.*B2_pattern(2,:),...
        1./barLength_square.*B2_pattern(3,:),...
        1./barLength_square.*B2_pattern(4,:),...
        1./barLength_square.*B2_pattern(5,:),...
        1./barLength_square.*B2_pattern(6,:)];
        
    U_temp=[U(NodeIndex1_temp,:) U(NodeIndex2_temp,:)];
    
    B2_U=[dot(B2_temp(:,1:6),U_temp,2),...
          dot(B2_temp(:,7:12),U_temp,2),...
          dot(B2_temp(:,13:18),U_temp,2),...
          dot(B2_temp(:,19:24),U_temp,2),...
          dot(B2_temp(:,25:30),U_temp,2),...
          dot(B2_temp(:,31:36),U_temp,2)];
      
    B1_B2_U = B1_temp + B2_U;
    
    K_temp_part1=[B1_B2_U(:,1).*B1_B2_U,...
                  B1_B2_U(:,2).*B1_B2_U,...
                  B1_B2_U(:,3).*B1_B2_U,...
                  B1_B2_U(:,4).*B1_B2_U,...
                  B1_B2_U(:,5).*B1_B2_U,...
                  B1_B2_U(:,6).*B1_B2_U];
              
    K_temp=C.*barArea.*barLength.*K_temp_part1+Sx.*barArea.*barLength.*B2_temp;    

    index1=3*NodeIndex1_temp-2;
    index2=3*NodeIndex2_temp-2;

    index_dim=[index1, index1+1, index1+2,...
         index2, index2+1, index2+2];
    
    index=zeros(N_bar,6,6,2);
    
    index(:,1,:,2)=index_dim;    
    index(:,2,:,2)=index_dim;
    index(:,3,:,2)=index_dim;
    index(:,4,:,2)=index_dim;
    index(:,5,:,2)=index_dim;
    index(:,6,:,2)=index_dim;
    
    index(:,:,1,1)=index_dim;    
    index(:,:,2,1)=index_dim;
    index(:,:,3,1)=index_dim;
    index(:,:,4,1)=index_dim;
    index(:,:,5,1)=index_dim;
    index(:,:,6,1)=index_dim;
    
    index=reshape(index, [N_bar*36,2]) ;
    
    K_temp=K_temp(:);
    length(K_temp);%Nb*36
    for i=1:length(K_temp)
        Kbar(index(i,1),index(i,2))=Kbar(index(i,1),index(i,2))+K_temp(i);
    end
    
end