

[Ru, Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
[Rd, Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
[RT, KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
[~, Kud] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'ud');
[~, Kdu] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');
% [~, KuT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
%     constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
[~, KdT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'dT');
[~, KTd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'Td');
Ru_active = Ru(active_indices_u);
Rd_active = Rd(active_indices_d);
RT_active = RT(active_indices_T);
Ku_active = Ku(active_indices_u, active_indices_u);
Kd_active = Kd(active_indices_d, active_indices_d);
KT_active = KT(active_indices_T, active_indices_T);
Kud_active = Kud(active_indices_u, active_indices_d);
Kdu_active = Kdu(active_indices_d, active_indices_u);
KdT_active = KdT(active_indices_d, active_indices_T);
KTd_active = KTd(active_indices_T, active_indices_d);
KuT_active = KuT(active_indices_u, active_indices_T);
delta_d = (Kd_active-KdT_active*(KT_active\KTd_active)-Kdu_active*(Ku_active\Kud_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\KTd_active))...
    \(Rd_active-Kdu_active*(Ku_active\Ru_active)-KdT_active*(KT_active\RT_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\RT_active));
delta_T = -KT_active\(RT_active+KTd_active*delta_d);
Sol_T(active_indices_T) = Sol_T(active_indices_T) + delta_T;
T_times = T_times + 1;

for uditer = 1:100
    Ru= Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
    Rd= Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
    Ru_active = Ru(active_indices_u);
    Rd_active = Rd(active_indices_d);
    if norm(Ru_active)<= Ru_tol * Ru_ref&&norm(Rd_active)<= Rd_tol * Rd_ref
        break
    end
    
    %½»Ìæ·¨
%     nt_u_converged = false;
%     for nt_u=1:max_iterT
%         [Ru,Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
%             constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
%         Ku_active = Ku(active_indices_u, active_indices_u);
%         Ru_active = Ru(active_indices_u);
%         fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
%         %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
%         if norm(Ru_active)<= R_tol * Ru_tol
%             Nt_converged = true;break;
%         end
%         Sol_u(active_indices_u) = Sol_u(active_indices_u) - Ku_active\Ru_active;
%         %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
%         u_times = u_times + 1;
%     end
%     if ~nt_u_converged
%         converged = false;
%     end
    
    %Êæ¶û²¹·¨
    [~, Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
    [~, Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');

    [~, Kud] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'ud');
    [~, Kdu] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');

    Ku_active = Ku(active_indices_u, active_indices_u);
    Kd_active = Kd(active_indices_d, active_indices_d);

    Kud_active = Kud(active_indices_u, active_indices_d);
    %         Kud_active = zeros(sum(active_indices_u), sum(active_indices_d));
    Kdu_active = Kdu(active_indices_d, active_indices_u);
    su = (Ku_active-Kud_active*(Kd_active\Kdu_active))\(Ru_active-Kud_active*(Kd_active\Rd_active));
%         sd = (Kd_active-KdT_active*(KT_active\KTd_active)-Kdu_active*(Ku_active\Kud_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\KTd_active))...
%             \(Rd_active-Kdu_active*(Ku_active\Ru_active)-KdT_active*(KT_active\RT_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\RT_active));
%         %åŒ–ç®€sd
%         %sd = (Kdu_active*Ku_active\(KuT_active*KT_active\RT_active-Kud_active)+)\();
%         su = -Ku_active\(Ru_active+(Kud_active-KuT_active*(KT_active\KTd_active))*sd-KuT_active*(KT_active\RT_active));
    
    nt_d_converged = false;
    for nt_d=1:max_iterT
        [Rd,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
        Kd_active = Kd(active_indices_d, active_indices_d);
        Rd_active = Rd(active_indices_d);
        fprintf('L2-norm of Rd = %g\n',norm(Rd_active));
        %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
        if norm(Rd_active)<= Rd_tol * Rd_ref
            Nt_converged = true;break;
        end
        Sol_d(active_indices_d) = Sol_d(active_indices_d) - Kd_active\Rd_active;
        %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
        d_times = d_times + 1;
    end
    if ~nt_d_converged
        converged = false;
    end
end