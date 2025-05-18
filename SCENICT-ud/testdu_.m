s_u_chushi =Sol_u;
s_u_raodong=+0.3*ones(length(s_u_chushi),1);
s_u_gai=s_u_chushi+s_u_raodong;
[~,Kdu] = Assemble_half_cohesive(s_u_gai,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');
deltard=Kdu*s_u_raodong;
[Rd1] = Assemble_half_cohesive(s_u_gai,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
Rd+deltard;