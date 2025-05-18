s_u_chushi =Sol_u(active_indices_u);%初始节点温度
s_u_raodong=-0.5*ones(length(s_u_chushi),1);%给节点温度扰动
s_u_gai = Sol_u;
s_u_gai(active_indices_u)=s_u_chushi+s_u_raodong;%改变后的节点温度
s_u_gai(~active_indices_u) = Sol_u(~active_indices_u);
[~,Ku] = Assemble_half_cohesive(s_u_gai,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
deltaru=Ku(active_indices_u,active_indices_u)*s_u_raodong;%
Ru_right = Assemble_half_cohesive(s_u_gai,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
Ru_right = Ru_right(active_indices_u);
Ru_test = Ru(active_indices_u)+deltaru;