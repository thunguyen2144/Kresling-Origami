function  [Mspr,Cspr]= Solve_Moment(obj,theta)

    theta1=obj.theta1;
    theta2=obj.theta2;

    sprRotK=obj.rot_spr_K_vec;
    spr_StressFree=obj.theta_stress_free_vec;

    Mspr=zeros(size(theta));
    Cspr=zeros(size(theta));

    A=size(theta);
    N_rotSpr=A(1);

    % This following part is not vectorized yet
    % The logical judgement for facade penertration is not easily 
    % vectorized due to the multiple "if" statements. However, 
    % because this part is not expensive computationally, we can
    % safely live with the for loop. 
    
    for i=1:N_rotSpr            
        if theta(i)<theta1
            ratio=(sec(pi*(theta(i)-theta1)/2/theta1))^2;
            Cspr(i)=sprRotK(i)*ratio;
            Mspr(i)=sprRotK(i)*(theta(i)-spr_StressFree(i))...
                +(2*theta1*sprRotK(i)/pi)*...
                tan(pi*(theta(i)-theta1)/2/theta1)-...
                sprRotK(i)*theta(i)+sprRotK(i)*theta1;
     
        elseif theta(i)>theta2
            ratio=(sec(pi*(theta(i)-theta2)/(4*pi-2*theta2)))^2;
            Cspr(i)=sprRotK(i)*ratio;
            Mspr(i)=sprRotK(i)*(theta(i)-spr_StressFree(i))...
                +(2*(2*pi-theta2)*sprRotK(i)/pi)...
                *tan(pi*(theta(i)-theta2)/(4*pi-2*theta2))-...
                sprRotK(i)*theta(i)+sprRotK(i)*theta2;
        else   
            Cspr(i)=sprRotK(i); 
            Mspr(i)=sprRotK(i)*(theta(i)-spr_StressFree(i));
        end

    end
        
end
