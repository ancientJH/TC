function [Sol_1,Sol_2,a_sctimes,a_nttimes,b_nttimes] = doublefieldschur(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf,Ru_ref,Rd_ref,RT_ref,...
        active_indices_u,active_indices_d,active_indices_T,Ru_tol,Rd_tol,RT_tol)

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
a_nttimes = 0;
a_sctimes = 0;
b_nttimes = 0;
for doubleiter=1:100
    
    R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
    R2 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
    R1_active = R1(active_indices_1);
    R2_active = R2(active_indices_2);
    
%     if doubleiter == 1
%         R1_ref = norm(R1_active);
%         R2_ref = norm(R2_active);
%     end
    fprintf('L2-norm of R1 = %g, L2-norm of R2 = %g\n',norm(R1_active),norm(R2_active));
    if norm(R1_active)<= R1_tol * R1_ref&&norm(R2_active)<= R2_tol * R2_ref
        break
    end
    if doubleiter == 1
        last_R1ref = norm(R1_active);
        R1_sc_ref = norm(R1_active);
        thisflag='sc法';
    else
        thisflag=nextflag;
    end
    if thisflag=='s求解'
        for nt_1=1:100
            [R1,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
            K1_active = K1(active_indices_1, active_indices_1);
            R1_active = R1(active_indices_1);
            fprintf('L2-norm of R1 = %g\n',norm(R1_active));
            if norm(R1_active)<= R1_tol * R1_ref
                nt1_converged = true;break;
            end
            Sol_1(active_indices_1) = Sol_1(active_indices_1) - K1_active\R1_active;
            a_nttimes = a_nttimes + 1;
            if whichhalf(1)=='u'
                Sol_u=Sol_1;
            elseif whichhalf(1)=='d'
                Sol_d=Sol_1;
            elseif whichhalf(1)=='T'
                Sol_T=Sol_1;
            end

            thisflag = 's求解';
            nextflag='sc法';
        end
    elseif thisflag=='sc法'

        [R1,K1] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
        [R2,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
        [~,K12] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf);
        [~,K21] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,reverse(whichhalf));
        %             [~,KuT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        %                 constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
        R1_active = R1(active_indices_1);
        R2_active = R2(active_indices_2);
        K1_active = K1(active_indices_1, active_indices_1);
        K2_active = K2(active_indices_2, active_indices_2);
        K12_active = K12(active_indices_1, active_indices_2);
        %             K21_active = zeros(sum(active_indices_2), sum(active_indices_1));
        K21_active = K21(active_indices_2, active_indices_1);
        %             KuT_active = KuT(active_indices_u, active_indices_T);
        for x= 1:-0.1:0
            try
                chol(K1_active-x^2*K12_active*(K2_active\K21_active));
                delta_1 = -(K1_active-x^2*K12_active*(K2_active\K21_active))\(R1_active-x*K12_active*(K2_active\R2_active));
                Sol_1(active_indices_1) = Sol_1(active_indices_1) + delta_1;
                a_sctimes = a_sctimes + 1;
                thisflag = 'sc法';
                if whichhalf(1)=='u'
                    Sol_u=Sol_1;
                elseif whichhalf(1)=='d'
                    Sol_d=Sol_1;
                elseif whichhalf(1)=='T'
                    Sol_T=Sol_1;
                end
                break
            catch
                continue
            end
        end


    end
%             delta_1 = -(K1_active-K12_active*(K2_active\K21_active))\(R1_active);

    nt2_converged = false;
    for nt_2=1:100
        [R2,K2] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(2));
        K2_active = K2(active_indices_2, active_indices_2);
        R2_active = R2(active_indices_2);
        fprintf('L2-norm of R2 = %g\n',norm(R2_active));
        if norm(R2_active)<= R2_tol * R2_ref
            nt2_converged = true;break;
        end
        Sol_2(active_indices_2) = Sol_2(active_indices_2) - K2_active\R2_active;
        b_nttimes = b_nttimes + 1;
        if whichhalf(2)=='u' 
            Sol_u=Sol_2;
        elseif whichhalf(2)=='d'
            Sol_d=Sol_2;
        elseif whichhalf(2)=='T'
            Sol_T=Sol_2;
        end
    end
    if ~nt2_converged
        converged = false;
        error('牛顿法失败')
    end
    
    if thisflag == 'sc法'
        R1 = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,whichhalf(1));
        R1_active = R1(active_indices_1);
        last_R1ref=R1_sc_ref;
        R1_sc_ref=norm(R1_active);

        if R1_sc_ref>last_R1ref
            %         if mod(a_sctimes+1,15)==0
            %             nextflag='s求解';
            nextflag='s求解';

        else
            nextflag='sc法';
        end
    end
end


