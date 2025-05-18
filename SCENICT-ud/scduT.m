

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
[~, KuT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
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
% delta_d = (Kd_active-KdT_active*(KT_active\KTd_active)-Kdu_active*(Ku_active\Kud_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\KTd_active))...
%     \(Rd_active-Kdu_active*(Ku_active\Ru_active)-KdT_active*(KT_active\RT_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\RT_active));
delta_d = (Kd_active-Kdu_active*(Ku_active\Kud_active))\(Rd_active-Kdu_active*(Ku_active\Ru_active));
Sol_d(active_indices_d) = Sol_d(active_indices_d) + delta_d;
d_times = d_times + 1;
for uTiter = 1:100
    Ru = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
    RT = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
    Ru_active = Ru(active_indices_u);
    
    RT_active = RT(active_indices_T);
    if norm(Ru_active)<= R_tol * Ru_tol&&norm(RT_active)<= R_tol * RT_tol
        break
    end
    nt_u_converged = false;
    for nt_u=1:max_iterT
        [Ru,Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
        Ku_active = Ku(active_indices_u, active_indices_u);
        Ru_active = Ru(active_indices_u);
        fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
        %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
        if norm(Ru_active)<= R_tol * Ru_tol
            Nt_converged = true;break;
        end
        Sol_u(active_indices_u) = Sol_u(active_indices_u) - Ku_active\Ru_active;
        %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
        u_times = u_times + 1;
    end
    if ~nt_u_converged
        converged = false;
    end
    
    nt_T_converged = false;
    for nt_T=1:max_iterT
        [RT,KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
        KT_active = KT(active_indices_T, active_indices_T);
        RT_active = RT(active_indices_T);
        fprintf('L2-norm of RT = %g\n',norm(RT_active));
        %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
        if norm(RT_active)<= R_tol * RT_tol
            Nt_converged = true;break;
        end
        Sol_T(active_indices_T) = Sol_T(active_indices_T) - KT_active\RT_active;
        %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
        T_times = T_times + 1;
    end
    if ~nt_T_converged
        converged = false;
    end
end