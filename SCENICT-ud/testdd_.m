s_d_chushi =Sol_d;
s_d_raodong=+0.6*ones(length(s_d_chushi),1);
s_d_gai=s_d_chushi+s_d_raodong;
[~,Kd] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
deltard=Kd*s_d_raodong;
[Rd1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
Rd+deltard;