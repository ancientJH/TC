s_d_chushi =Sol_d;
s_d_raodong=+0.3*ones(length(s_d_chushi),1);
s_d_gai=s_d_chushi+s_d_raodong;
[Rd,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
[Rd1,Kd1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
[Ru,Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
[Ru1,Ku1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
[~,Kud] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'ud');
[~,Kdu] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');
[~,Kud1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'ud');
[~,Kdu1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');
% deltaru=Kud1*s_d_gai-Kud*s_d_chushi;
deltaru=Kud1*s_d_raodong;
[Ru1] = Assemble_half_cohesive(Sol_u,s_d_gai,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
Ru+deltaru;
% du=(Ku - Kud*(Kd\Kdu))\(Ru-Kud*(Kd\Rd));
% du1=(Ku1 - Kud1*(Kd1\Kdu1))\(Ru1-Kud1*(Kd1\Rd1));
% Rud=Ru-Ku*du;
% Rud1=Ru1-Ku1*du1;
% aaa=0.7*Rud+0.7*Kud*0.3*ones(2159,1);