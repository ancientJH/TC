s_t_chushi =Sol_T;
s_t_raodong=-100*ones(length(s_t_chushi),1);
s_t_gai=s_t_chushi+s_t_raodong;
[~,KTT] = Assemble_half_cohesive(Sol_u,Sol_d,s_t_gai,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
deltarT=KTT*s_t_raodong;
[RT1] = Assemble_half_cohesive(Sol_u,Sol_d,s_t_gai,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
RT+deltarT;