for iter_noall=1:max_iterT
        %for i = 1:2
        flag=0;
        Ru = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
        Rd = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
        RT = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
        Ru_active = Ru(active_indices_u);
        Rd_active = Rd(active_indices_d);
        RT_active = RT(active_indices_T);

        fprintf('L2-norm of Ru = %g, L2-norm of Rd = %g, L2-norm of RT = %g\n',norm(Ru_active),norm(Rd_active),norm(RT_active));
%         if norm(Rd_active) <= R_tol * sqrt(length(Rd_active)) && norm(RT_active) <= R_tol * sqrt(length(RT_active))/102.6
%             flag=1;
%         end

        if norm(Ru_active) <= R_tol * Ru_tol && norm(Rd_active) <= R_tol * Rd_tol && norm(RT_active) <= R_tol * RT_tol
            converged = true;
            break;
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



%dT迭代开始
for dt_iter = 1:max_iterT
    Rd = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
    RT = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
    Rd_active = Rd(active_indices_d);
    RT_active = RT(active_indices_T);
    fprintf('L2-norm of Rd = %g, L2-norm of RT = %g\n',norm(Rd_active),norm(RT_active));
    if norm(Rd_active) <= R_tol * Rd_tol && norm(RT_active) <= R_tol * RT_tol
        break;
    end

   nt_d_converged = false;
    for nt_d=1:max_iterT
        [Rd,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
        Kd_active = Kd(active_indices_d, active_indices_d);
        Rd_active = Rd(active_indices_d);
        fprintf('L2-norm of Rd = %g\n',norm(Rd_active));
        %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
        if norm(Rd_active)<= R_tol * Rd_tol
            Nt_converged = true;break;
        end
        Sol_d(active_indices_d) = Sol_d(active_indices_d) - Kd_active\Rd_active;
        %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
        d_times = d_times + 1;
    end
    if ~nt_d_converged
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

end