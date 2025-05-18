load_pre = load_pre + deltau;
Sol_T(BCIndices_T) = BCVal_T(BCIndices_T);
Sol_T_record = Sol_T;
switch iterationType
    case 'monolithic'
        Sol(BCIndices) = BCVal(BCIndices) * load_pre;
        %load steps are to ease the solution of nonlinear problems.
        % It can be taken as 1 or [0,1], or something like 0:0.1:1
        % The body force and traction are scaled by the same factor
    case 'staggering'
        Sol_u(BCIndices_u) = BCVal_u(BCIndices_u);
        % Sol_d(BCIndices_d) = BCVal_d(BCIndices_d) * load_steps(step_no);
        Sol_d(BCIndices_d) = BCVal_d(BCIndices_d);
end
for nt_T=1:max_iterT
    [RT,KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
    KT_active = KT(active_indices_T, active_indices_T);
    RT_active = RT(active_indices_T);
    fprintf('L2-norm of RT = %g\n',norm(RT_active));
    %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
    if norm(RT_active)<= 1e-6 * RT_ref
        Nt_converged = true;break;
    end
    Sol_T(active_indices_T) = Sol_T(active_indices_T) - KT_active\RT_active;
    %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
    ntT_times = ntT_times + 1;
end

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
[~, Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
[~, Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
[~, KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
Ku_active = Ku(active_indices_u, active_indices_u);
Kd_active = Kd(active_indices_d, active_indices_d);
KT_active = KT(active_indices_T, active_indices_T); 
Kud_active = Kud(active_indices_u, active_indices_d);
Kdu_active = Kdu(active_indices_d, active_indices_u);
KdT_active = KdT(active_indices_d, active_indices_T);
KTd_active = KTd(active_indices_T, active_indices_d);
KuT_active = KuT(active_indices_u, active_indices_T);
KTu_active = zeros(sum(active_indices_T), sum(active_indices_u));
uT = abs(sqrt(max(eigs((KT_active\KTu_active)*(Ku_active\KuT_active)))))
ud = abs(sqrt(max(eigs((Kd_active\Kdu_active)*(Ku_active\Kud_active)))))
Td = abs(sqrt(max(eigs((Kd_active\KdT_active)*(KT_active\KTd_active)))))