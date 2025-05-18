s_t_chushi =Sol_T;
s_t_raodong=-1*ones(length(s_t_chushi),1);
s_t_gai=s_t_chushi+s_t_raodong;
[~,KuT1] = Assemble_half_cohesive(Sol_u,Sol_d,s_t_gai,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
[~,KuT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
deltaru=KuT*s_t_raodong;
[Ru1] = Assemble_half_cohesive(Sol_u,Sol_d,s_t_gai,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
Ru+deltaru;
err = norm(ans-Ru1)/norm(Ru)