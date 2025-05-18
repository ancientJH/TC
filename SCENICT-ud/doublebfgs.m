function [Sol_1,Sol_2,a_bftimes,a_nttimes,b_nttimes] = doublebfgs(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf,Ru_ref,Rd_ref,RT_ref,...
    active_indices_u,active_indices_d,active_indices_T,Ru_tol,Rd_tol,RT_tol)

a_nttimes = 0;
a_bftimes = 0;
b_nttimes = 0;
if whichhalf(1)=='u'
    R1_tol = Ru_tol;
    R1_ref = Ru_ref;
    active_indices_1 = active_indices_u;
    Sol_1=Sol_u;
elseif whichhalf(1)=='d'
    R1_tol = Rd_tol;
    R1_ref = Rd_ref;
    active_indices_1 = active_indices_d;
    Sol_1=Sol_d;
elseif whichhalf(1)=='T'
    R1_tol = RT_tol;
    R1_ref = RT_ref;
    active_indices_1 = active_indices_T;
    Sol_1=Sol_T;
end
if whichhalf(2)=='u'
    R2_tol = Ru_tol;
    R2_ref = Ru_ref;
    active_indices_2 = active_indices_u;
    Sol_2=Sol_u;
elseif whichhalf(2)=='d'
    R2_tol = Rd_tol;
    R2_ref = Rd_ref;
    active_indices_2 = active_indices_d;
    Sol_2=Sol_d;
elseif whichhalf(2)=='T'
    R2_tol = RT_tol;
    R2_ref = RT_ref;
    active_indices_2 = active_indices_T;
    Sol_2=Sol_T;
end
Sol_active = zeros(sum(active_indices_1)+sum(active_indices_2),2);
R_active = zeros(sum(active_indices_1)+sum(active_indices_2),2);
Hs = zeros(sum(active_indices_1)+sum(active_indices_2),2*(sum(active_indices_2)+sum(active_indices_1)));
[~,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
[~,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
K1_active = K1(active_indices_1, active_indices_1);
K2_active = K2(active_indices_2, active_indices_2);
K_mono_inverse = [inv(K1_active),zeros(sum(active_indices_1),sum(active_indices_2));zeros(sum(active_indices_2),sum(active_indices_1)),inv(K2_active)];
Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2)) = K_mono_inverse;
Sol_active(:,1) = [Sol_1(active_indices_1);Sol_2(active_indices_2)];

R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
R2 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
R1_active = R1(active_indices_1);
R2_active = R2(active_indices_2);
R_active(:,1)=[R1_active;R2_active];
Sol_active(:,2) = Sol_active(:,1) - Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2))*R_active(:,1);
a_nttimes = a_nttimes + 1;
b_nttimes = b_nttimes + 1;
Sol_1(active_indices_1) = Sol_active(1:sum(active_indices_1),2);
if whichhalf(1)=='u'
    Sol_u=Sol_1;
elseif whichhalf(1)=='d'
    Sol_d=Sol_1;
elseif whichhalf(1)=='T'
    Sol_T=Sol_1;
end
Sol_2(active_indices_2) = Sol_active(sum(active_indices_1)+1:end,2);
if whichhalf(2)=='u'
    Sol_u=Sol_2;
elseif whichhalf(2)=='d'
    Sol_d=Sol_2;
elseif whichhalf(2)=='T'
    Sol_T=Sol_2;
end
R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
R2 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
R1_active = R1(active_indices_1);
R2_active = R2(active_indices_2);
R_active(:,2)=[R1_active;R2_active];


for doubleiter=1:100

    R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
    R2 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
    R1_active = R1(active_indices_1);
    R2_active = R2(active_indices_2);

    fprintf('L2-norm of R1 = %g, L2-norm of R2 = %g\n',norm(R1_active),norm(R2_active));
    if norm(R1_active)<= R1_tol * R1_ref&&norm(R2_active)<= R2_tol * R2_ref
        break
    end
    if mod(doubleiter,2)==1
        if mod(doubleiter,5)==0
            [R1,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            [R2,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
            K1_active = K1(active_indices_1, active_indices_1);
            K2_active = K2(active_indices_2, active_indices_2);
            R1_active = R1(active_indices_1);
            R2_active = R2(active_indices_2);
            Sol_active(:,1) = Sol_active(:,2) - [K1_active,zeros(sum(active_indices_1),sum(active_indices_2));zeros(sum(active_indices_2),sum(active_indices_1)),K2_active]\[R1_active;R2_active];
            a_nttimes = a_nttimes + 1;
            b_nttimes = b_nttimes + 1;
            Sol_1(active_indices_1) = Sol_active(1:sum(active_indices_1),1);
            if whichhalf(1)=='u'
                Sol_u=Sol_1;
            elseif whichhalf(1)=='d'
                Sol_d=Sol_1;
            elseif whichhalf(1)=='T'
                Sol_T=Sol_1;
            end
            Sol_2(active_indices_2) = Sol_active(sum(active_indices_1)+1:end,1);
            if whichhalf(2)=='u'
                Sol_u=Sol_2;
            elseif whichhalf(2)=='d'
                Sol_d=Sol_2;
            elseif whichhalf(2)=='T'
                Sol_T=Sol_2;
            end
            [R1,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            [R2,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
            K1_active = K1(active_indices_1, active_indices_1);
            K2_active = K2(active_indices_2, active_indices_2);
            R1_active = R1(active_indices_1);
            R2_active = R2(active_indices_2);
            R_active(:,1)=[R1_active;R2_active];
            K_mono_inverse = [inv(K1_active),zeros(sum(active_indices_1), sum(active_indices_2));zeros(sum(active_indices_2), sum(active_indices_1)),inv(K2_active)];
            Hs(1:size(K_mono_inverse,1),size(K_mono_inverse,2)+1:end) = K_mono_inverse;
 
        else
            delta_Sol = Sol_active(:,2) - Sol_active(:,1);
            delta_R = R_active(:,2) - R_active(:,1);
            Hs(1:size(K_mono_inverse,1),size(K_mono_inverse,2)+1:end) = (eye(size(K_mono_inverse,1)) - (delta_Sol*delta_R')/(delta_Sol'*delta_R))*Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2))...
                *(eye(size(K_mono_inverse,1)) - (delta_Sol*delta_R')/(delta_Sol'*delta_R))' + (delta_Sol*delta_Sol')/(delta_Sol'*delta_R);
            Sol_active(:,1) = Sol_active(:,2) - Hs(1:size(K_mono_inverse,1),size(K_mono_inverse,2)+1:end)*R_active(:,2);
            a_nttimes = a_nttimes + 1;
            b_nttimes = b_nttimes + 1;
            Sol_1(active_indices_1) = Sol_active(1:sum(active_indices_1),1);
            if whichhalf(1)=='u'
                Sol_u=Sol_1;
            elseif whichhalf(1)=='d'
                Sol_d=Sol_1;
            elseif whichhalf(1)=='T'
                Sol_T=Sol_1;
            end
            Sol_2(active_indices_2) = Sol_active(sum(active_indices_1)+1:end,1);
            if whichhalf(2)=='u'
                Sol_u=Sol_2;
            elseif whichhalf(2)=='d'
                Sol_d=Sol_2;
            elseif whichhalf(2)=='T'
                Sol_T=Sol_2;
            end
            R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            R2 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
            R1_active = R1(active_indices_1);
            R2_active = R2(active_indices_2);
            R_active(:,1)=[R1_active;R2_active];
        end
    else
        if mod(doubleiter,5)==0
            [R1,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            [R2,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
            K1_active = K1(active_indices_1, active_indices_1);
            K2_active = K2(active_indices_2, active_indices_2);
            R1_active = R1(active_indices_1);
            R2_active = R2(active_indices_2);
            Sol_active(:,2) = Sol_active(:,1) - [K1_active,zeros(sum(active_indices_1), sum(active_indices_2));zeros(sum(active_indices_2), sum(active_indices_1)),K2_active]\[R1_active;R2_active];
            a_nttimes = a_nttimes + 1;
            b_nttimes = b_nttimes + 1;
            Sol_1(active_indices_1) = Sol_active(1:sum(active_indices_1),2);
            if whichhalf(1)=='u'
                Sol_u=Sol_1;
            elseif whichhalf(1)=='d'
                Sol_d=Sol_1;
            elseif whichhalf(1)=='T'
                Sol_T=Sol_1;
            end
            Sol_2(active_indices_2) = Sol_active(sum(active_indices_1)+1:end,2);
            if whichhalf(2)=='u'
                Sol_u=Sol_2;
            elseif whichhalf(2)=='d'
                Sol_d=Sol_2;
            elseif whichhalf(2)=='T'
                Sol_T=Sol_2;
            end
            [R1,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            [R2,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
            K1_active = K1(active_indices_1, active_indices_1);
            K2_active = K2(active_indices_2, active_indices_2);
            R1_active = R1(active_indices_1);
            R2_active = R2(active_indices_2);
            R_active(:,2)=[R1_active;R2_active];
            K_mono_inverse = [inv(K1_active),zeros(sum(active_indices_1), sum(active_indices_2));zeros(sum(active_indices_2), sum(active_indices_1)),inv(K2_active)];
            Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2)) = K_mono_inverse;

        else
            delta_Sol = Sol_active(:,1) - Sol_active(:,2);
            delta_R = R_active(:,1) - R_active(:,2);
            Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2)) = (eye(size(K_mono_inverse,1)) - (delta_Sol*delta_R')/(delta_Sol'*delta_R))*Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2))...
                *(eye(size(K_mono_inverse,1)) - (delta_Sol*delta_R')/(delta_Sol'*delta_R))' + (delta_Sol*delta_Sol')/(delta_Sol'*delta_R);
            Sol_active(:,2) = Sol_active(:,1) - Hs(1:size(K_mono_inverse,1),1:size(K_mono_inverse,2))*R_active(:,1);
            a_nttimes = a_nttimes + 1;
            b_nttimes = b_nttimes + 1;
            Sol_1(active_indices_1) = Sol_active(1:sum(active_indices_1),2);
            if whichhalf(1)=='u'
                Sol_u=Sol_1;
            elseif whichhalf(1)=='d'
                Sol_d=Sol_1;
            elseif whichhalf(1)=='T'
                Sol_T=Sol_1;
            end
            Sol_2(active_indices_2) = Sol_active(sum(active_indices_1)+1:end,2);
            if whichhalf(2)=='u'
                Sol_u=Sol_2;
            elseif whichhalf(2)=='d'
                Sol_d=Sol_2;
            elseif whichhalf(2)=='T'
                Sol_T=Sol_2;
            end
            R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            R2 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
            R1_active = R1(active_indices_1);
            R2_active = R2(active_indices_2);
            R_active(:,2)=[R1_active;R2_active];

        end
    end
    
end

