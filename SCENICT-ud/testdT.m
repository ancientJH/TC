s_T_chushi =Sol_T;
s_T_raodong=1000*ones(length(s_T_chushi),1);
s_T_gai=s_T_chushi+s_T_raodong;
[~,KdT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'dT');
[~,KdT1] = Assemble_half_cohesive(Sol_u,Sol_d,s_T_gai,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'dT');
deltard=KdT1*s_T_raodong;
[Rd1] = Assemble_half_cohesive(Sol_u,Sol_d,s_T_gai,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
Rd+deltard;
err = norm(ans-Rd1)/norm(Rd)