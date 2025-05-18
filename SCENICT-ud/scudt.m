for iter_noSc=1:max_iterT
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


    if norm(Ru_active) >= Ru_tol||mod(u_times+1,4)==0
        % if norm(Ru_active)>=R_tol * sqrt(length(Ru_active))/(0.0021*5e-3)||mod(u_times,4)==0
        [~, Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
        Ku_active = Ku(active_indices_u, active_indices_u);
        su = -Ku_active\Ru_active;
        Sol_u(active_indices_u) = Sol_u(active_indices_u) + su;
        u_times = u_times + 1;
    else
        % a=(eye(length(Ku_active))+(Ku_active\Kud_active)*inv((eye(length(Kd_active))-(Kd_active\KdT_active)*(KT_active\KTd_active)))*(Kd_active\Kdu_active)...
        %     +(Ku_active\KuT_active)*(KT_active\KTd_active)*inv((eye(length(Kd_active))-(Kd_active\KdT_active)*(KT_active\KTd_active)))*(Kd_active\Kdu_active));
        % b=(-Ku_active*Ru_active-(Ku_active\Kud_active)*inv((eye(length(Kd_active))-(Kd_active\KdT_active)*(KT_active\KTd_active)))*((Kd_active\KdT_active)*(KT_active\RT_active)-(Kd_active\Rd_active))...
        %     +(Ku_active\KuT_active)*(KT_active\(RT_active+KTd_active*((eye(length(Kd_active))-(Kd_active\KdT_active)*(KT_active\KTd_active))\((Kd_active\KdT_active)*(KT_active\RT_active)-Kd_active\Rd_active)))));
        % su = a\b;
        % su = lsqminnorm((eye(length(Ku_active))+Ku_active\(KuT_active*(KT_active\KTd_active)...
        %     -Kud_active)*((Kd_active-KdT_active*(KT_active\KTd_active))\Kdu_active))...
        %     ,(-Ku_active\Ru_active+Ku_active\(KuT_active*(KT_active\KTd_active)-Kud_active)...
        %     *((Kd_active-KdT_active*(KT_active\KTd_active))\(-Rd_active+KdT_active*(KT_active\RT_active)))...
        %     +(Ku_active\KuT_active)*(KT_active\RT_active)));

        [~, Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
        [~, Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
        [~, KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
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
        Ku_active = Ku(active_indices_u, active_indices_u);
        Kd_active = Kd(active_indices_d, active_indices_d);
        KT_active = KT(active_indices_T, active_indices_T);
        Kud_active = Kud(active_indices_u, active_indices_d);
        %         Kud_active = zeros(sum(active_indices_u), sum(active_indices_d));
        Kdu_active = Kdu(active_indices_d, active_indices_u);
        %         Kdu_active = zeros(sum(active_indices_d), sum(active_indices_u));
        KuT_active = KuT(active_indices_u, active_indices_T);
        %         KuT_active = zeros(sum(active_indices_u), sum(active_indices_T));
        KdT_active = KdT(active_indices_d, active_indices_T);
        %         KdT_active = zeros(sum(active_indices_d), sum(active_indices_T));
        KTd_active = KTd(active_indices_T, active_indices_d);
        %         KTd_active = zeros(sum(active_indices_T), sum(active_indices_d));
        sd = (Kd_active-KdT_active*(KT_active\KTd_active)-Kdu_active*(Ku_active\Kud_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\KTd_active))...
            \(Rd_active-Kdu_active*(Ku_active\Ru_active)-KdT_active*(KT_active\RT_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\RT_active));
        %化简sd
        %sd = (Kdu_active*Ku_active\(KuT_active*KT_active\RT_active-Kud_active)+)\();
        su = -Ku_active\(Ru_active+(Kud_active-KuT_active*(KT_active\KTd_active))*sd-KuT_active*(KT_active\RT_active));
        %         sT = -KT_active\(RT_active+KTd_active*sd);
        Sol_u(active_indices_u) = Sol_u(active_indices_u) + su;
        u_times = u_times + 1;
        %         Sol_d(active_indices_d) = Sol_d(active_indices_d) + sd;
        %         d_times = d_times + 1;
        %         Sol_T(active_indices_T) = Sol_T(active_indices_T) + sT;
        %         T_times = T_times + 1;

        %         su = (eye(length(Ku_active))+Ku_active\(KuT_active*(KT_active\KTd_active)...
        %             -Kud_active)*((Kd_active-KdT_active*(KT_active\KTd_active))\Kdu_active))...
        %             \(-Ku_active\Ru_active+Ku_active\(KuT_active*(KT_active\KTd_active)-Kud_active)...
        %             *((Kd_active-KdT_active*(KT_active\KTd_active))\(-Rd_active+KdT_active*(KT_active\RT_active)))...
        %             +(Ku_active\KuT_active)*(KT_active\RT_active));

    end
    %         su = inv(eye(length(Ku_active))+inv(Ku_active)*(KuT_active*inv(KT_active)*KTd_active-Kud_active)...
    %             *inv(Kd_active-KdT_active*inv(KT_active)*KTd_active)*Kdu_active)...
    %             *(-inv(Ku_active)*Ru_active+inv(Ku_active)*(KuT_active*inv(KT_active)*KTd_active-Kud_active)...
    %             *inv(Kd_active-KdT_active*inv(KT_active)*KTd_active)*(-Rd_active+KdT_active*inv(KT_active)*RT_active)...
    %             +inv(Ku_active)*KuT_active*inv(KT_active)*RT_active);
    %         su = -Ku_active\Ru_active;

    %         sd = (Kd_active-KdT_active*(KT_active\KTd_active)-Kdu_active*(Ku_active\Kud_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\KTd_active))...
    %             \(Rd_active-Kdu_active*(Ku_active\Ru_active)-KdT_active*(KT_active\RT_active)+Kdu_active*(Ku_active\KuT_active)*(KT_active\RT_active));
    %         su = -Ku_active*(Ru_active+(Kud_active-KuT_active*(KT_active\KTd_active))*sd-KuT_active*(KT_active\RT_active));

    %     end


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
        if norm(Rd_active)>Rd_tol||mod(d_times+1,4) == 0
            [~,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
            Kd_active = Kd(active_indices_d, active_indices_d);

            sd = -Kd_active\Rd_active;
            Sol_d(active_indices_d) = Sol_d(active_indices_d) + sd;
            
            d_times = d_times + 1;

        else
            [~, KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
            [~, KdT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'dT');
            [~, KTd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'Td');
            [~,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
            Kd_active = Kd(active_indices_d, active_indices_d);
            KT_active = KT(active_indices_T, active_indices_T);
            KdT_active = KdT(active_indices_d, active_indices_T);
            KTd_active = KTd(active_indices_T, active_indices_d);
            sd = -(Kd_active-KdT_active*(KT_active\KTd_active))\(Rd_active-KdT_active*(KT_active\RT_active));
            Sol_d(active_indices_d) = Sol_d(active_indices_d) + sd;
            h = h + 1;
            d_times = d_times + 1;
            %    d_times = d_times + 1;%Increment the number of displacement field iterations by one.
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