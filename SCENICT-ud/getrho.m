% 假设 load_steps 是一个数组，包含了您想要转换的小数
folderPath = 'D:\jjh\i7三场\_20241216ell0-0075\非绝对值容差\绝对值容差SCENICT-ud\Thermo8467ele';
loads = 0.0001:0.0001:0.01;

% 遍历 load_steps 数组，生成每个小数对应的文件名
for step_noo = 1:length(loads)
    % 获取当前步骤的小数值
    step_value = loads(step_noo);
    step_str = strrep(num2str(step_value, '%.5f'), '.', '_');
    filename = ['PF-', step_str, '.mat'];
    filePath = fullfile(folderPath, filename);
    load(filePath);

    if step_noo == 1
        Psi_plus_old = zeros(size(D_elem,1),size(D_elem,2));
        Psi_plus_rec = zeros(size(D_elem,1),size(D_elem,2));
    else
        step_value_last = loads(step_noo-1);
        step_str_last = strrep(num2str(step_value_last, '%.5f'), '.', '_');
        filename_last = ['PF-', step_str_last, '.mat'];
        filePath_last = fullfile(folderPath, filename_last);
        load(filePath_last,'Psi_plus_old','Psi_plus_rec');
    end
    % 显示生成的文件名
    
    disp(filePath);
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
    [~,KuT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'uT');
    [~,KTd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'Td');
    [~,KdT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'dT');
    Ku_active = Ku(active_indices_u, active_indices_u);
    Kd_active = Kd(active_indices_d, active_indices_d);
    Kud_active = Kud(active_indices_u, active_indices_d);
    Kdu_active = Kdu(active_indices_d, active_indices_u);
    KdT_active = KdT(active_indices_d, active_indices_T);
    KTd_active = KTd(active_indices_T, active_indices_d);
    KuT_active = KuT(active_indices_u, active_indices_T);
    KT_active = KT(active_indices_T, active_indices_T);
    KTu_active = sparse(sum(active_indices_T),sum(active_indices_u));

    K11 = KT_active; K12 = [KTu_active,KTd_active]; K21 = [KuT_active;KdT_active]; K22 = [Ku_active,Kud_active;Kdu_active,Kd_active];

    rhoa(step_noo,1) = abs(eigs([sparse(size(K11,1),size(K11,2)),K12; ...
        sparse(size(K21,1),size(K21,2)),sparse(size(K22,1),size(K22,2))], ...
        [K11,sparse(size(K12,1),size(K12,2));K21,K22],1));

    K11 = Ku_active; K12 = Kud_active; K21 = Kdu_active; K22 = Kd_active;

    rhob(step_noo,1) = abs(eigs([sparse(size(K11,1),size(K11,2)),K12; ...
        sparse(size(K21,1),size(K21,2)),sparse(size(K22,1),size(K22,2))], ...
        [K11,sparse(size(K12,1),size(K12,2));K21,K22],1));

end
save rho.mat