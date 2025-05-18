[~,Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
[~,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
[~,KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
[~,Kud] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'ud');
[~,Kdu] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');
[~,KdT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'dT');
[~,KTd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'Td');
[~,KuT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
% [~,KTu] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
%     constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'Tu');
Ku_active = Ku(active_indices_u, active_indices_u);
Kd_active = Kd(active_indices_d, active_indices_d);
KT_active = KT(active_indices_T, active_indices_T);
Kud_active = Kud(active_indices_u, active_indices_d);
Kdu_active = Kdu(active_indices_d, active_indices_u);
KdT_active = KdT(active_indices_d, active_indices_T);
KTd_active = KTd(active_indices_T, active_indices_d);
KuT_active = KuT(active_indices_u, active_indices_T);
KTu_active = zeros(sum(active_indices_T), sum(active_indices_u));



rhoT_ud = abs(eigs([KTu_active, KTd_active]*([Ku_active,Kud_active;Kdu_active,Kd_active]\[KuT_active;KdT_active]),KT_active,1))
rhou_Td = abs(eigs([KuT_active, Kud_active]*([KT_active,KTd_active;KdT_active,Kd_active]\[KTu_active;Kdu_active]),Ku_active,1));
rhod_Tu = abs(eigs([KdT_active, Kdu_active]*([KT_active,KTu_active;KuT_active,Ku_active]\[KTd_active;Kud_active]),Kd_active,1));


rhoT_ud = max(abs(eigs(KT_active\[KTu_active,KTd_active]*([Ku_active,Kud_active;Kdu_active,Kd_active]\[KuT_active;KdT_active]))));
rhou_Td = abs(eigs([KuT_active, Kud_active]*([KT_active,KTd_active;KdT_active,Kd_active]\[KTu_active;Kdu_active]),Ku_active,1));

rho_ud = abs(eigs(Kud_active*(Kd_active\Kdu_active),Ku_active,1))
rho_Td = abs(eigs(KTd_active*(Kd_active\KdT_active),KT_active,1))

%把u放在外面，Td之间交替法