s_d_chushi =Sol_d;
s_d_raodong=+0.05*ones(length(s_d_chushi),1);
s_d_gai=s_d_chushi+s_d_raodong;
[~,KTd1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'Td');
deltarT=KTd1*s_d_raodong;
[RT1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
RT+deltarT;
err = norm(ans-RT1)/norm(RT)