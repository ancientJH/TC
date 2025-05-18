%%
close all;
clear all;
clc;
if(exist('Thermo8467ele','dir') ~= 0)
    error ('Thermo8467ele/*.*');
else
    mkdir('Thermo8467ele');
end
elementType = 'P12D';
ElementType = 'P12D';% 'P12D' means 2D triangular element,'P22D' means 2D triangular element with 2-order shape function
problemType = 'temperature';         % Other options: 'Elasticity''PhaseField'

[nDim, nDoF, nNodes, nElements, nNodesElement, Coord, nEquations, ...
    DirichletBCs,IEN, LM, constitutive, body_force, traction] = ...
    ProblemDefinition(elementType, problemType);

% mesh
%  figure(1)
% X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
% Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
% %D_elem = ones(size(X));
% D_elem = linspace(1,size(IEN,2),size(IEN,2));
% 
% % Patch
% hold on
% XYZ=patch(X,Y,D_elem);
% axis equal;
% axis off;
% colorbar;


%generate gauss points on each element
% gauss_on_ele_y = zeros(size(IEN));
% gauss_on_ele_x = zeros(size(IEN));
switch ElementType
    case 'P12D'
        nQuad = 3; nQuadBdry = 2;   %Nquad is the nodes'number of one element,nQuadRdry is the nodes'number of one side
    case 'Q12D'
        nQuad = 4; nQuadBdry = 2;
end

% --------------------------------------------------------------------------
% Phase 2: Solve with load stepping
% --------------------------------------------------------------------------
% load_steps = -load_steps;
% figure(2)
% plot(load_steps);
% hold on;
% plot(cyc_record);
% title(sprintf('The load steps'));
% xlabel('time step');
% ylabel('load value');


global deltau T0 
T0 = 573;    %inital temperature

Sol = zeros(nEquations, 1);
[n,bc_check,BCIndices, BCVal,gc_vec] = DirichletBCs(Coord);   %boundary condition


if strcmp(problemType, 'Poisson') && sum(BCIndices) == 0 % fix the value of any node in case of no Dirichlet BC
    BCIndices(1) = true;
    BCVal(1) = 0;
end
active_indices = ~BCIndices;   %Not on the Dirichlet boundary, i.e., free boundary
% R_tol = 1e-14;    %----------->absolute tolerance of the residual
% R_tol = 1e-40; 
iterationType = 'staggering'; % Possible values: 'monolithic'or 'staggering'

if strcmp(iterationType, 'staggering')
    if strcmp(problemType, 'Elasticity') || strcmp(problemType, 'Poisson')
        iterationType = 'monolithic';
    else
        if strcmp(problemType,'temperature')
            d_indices = (nDim+1):(nDim+2):length(Sol); % Global indices in Sol
            T_indices = (nDim+2):(nDim+2):length(Sol);
            u_indices = setdiff(1:length(Sol), d_indices);
            u_indices = setdiff(u_indices,T_indices);
            BCIndices_u = BCIndices(u_indices); BCVal_u = BCVal(u_indices); % Indices of the u-subproblem
            BCIndices_d = BCIndices(d_indices); BCVal_d = BCVal(d_indices); % Indices of the d-subproblem
            BCIndices_T = BCIndices(T_indices); BCVal_T = BCVal(T_indices);
            active_indices_u = ~BCIndices_u; % Indices of the u-subproblem
            active_indices_d = ~BCIndices_d; % Indices of the d-subproblem
            active_indices_T = ~BCIndices_T;
            Sol(T_indices) = T0;
            Sol_u = Sol(u_indices); %initial solution of u
            Sol_d = Sol(d_indices); %initial solution of d
            Sol_T = Sol(T_indices); %initial solution of T
            LM_u = reshape(repmat(reshape((IEN-1)*nDim, 1, []), nDim, 1) ...
                + repmat((1:nDim)', 1, numel(IEN)), nNodesElement * nDim, []);
            LM_d = IEN;
            LM_T = IEN;

        end

        if strcmp(problemType,'phasefield')
            d_indices = (nDim+1):(nDim+1):length(Sol); % Global indices in Sol
            u_indices = setdiff(1:length(Sol), d_indices); % Global indices in Sol
            BCIndices_u = BCIndices(u_indices); BCVal_u = BCVal(u_indices); % Indices of the u-subproblem
            BCIndices_d = BCIndices(d_indices); BCVal_d = BCVal(d_indices); % Indices of the d-subproblem
            active_indices_u = ~BCIndices_u; % Indices of the u-subproblem
            active_indices_d = ~BCIndices_d; % Indices of the d-subproblem
            Sol_u = Sol(u_indices);
            Sol_d = Sol(d_indices);
            LM_u = reshape(repmat(reshape((IEN-1)*nDim, 1, []), nDim, 1) ...
                + repmat((1:nDim)', 1, numel(IEN)), nNodesElement * nDim, []);
            LM_d = IEN;
        end
    end
end

% ll=0;
 maxloadstep =50;
 pic_num0 = 1;     %draw dynamic picture


 global Psi_plus_rec Psi_plus_old %Psi_plus_old_d Psi_plus_rec_d Pure_Psi_plus driving_force
Psi_plus_old = zeros(nQuad,size(IEN,2));        %The maximum value of the Gaussian point psi+as of the last loading step  %nQuads is the number of gauss points.size(IEN,2) is the size of the elements.
Psi_plus_rec = zeros(nQuad,size(IEN,2));        % Maximum psi+
% Psi_plus_old_d = zeros(nQuad,size(IEN,2));      % g(d)psi+ of last step
% Psi_plus_rec_d = zeros(nQuad,size(IEN,2));      % g(d)psi+ of current step
% Pure_Psi_plus = zeros(nQuad,size(IEN,2));       % psi+ of current step
% driving_force = zeros(nQuad,size(IEN,2));       % driving force of phase field   g(d)*Psi_plus_rec


% elastic stress of gauss point
% global stressvector_gauss_y stressvector_gauss_x stressvector_gauss_xy Hydrostatic_stress 
% stressvector_gauss_y = zeros(nQuad,size(IEN,2));%nQuads is the number of gauss points.It is the y direction of the gauss points, you can change the direction by
% stressvector_gauss_x = zeros(nQuad,size(IEN,2));   
% stressvector_gauss_xy = zeros(nQuad,size(IEN,2));  
% Hydrostatic_stress = zeros(nQuad,size(IEN,2));
% 
% % Coord of gauss point
% global R_gauss x_gauss y_gauss strainvector_gauss_x strainvector_gauss_y strainvector_gauss_xy
% R_gauss = zeros(nQuad,size(IEN,2));
% x_gauss = zeros(nQuad,size(IEN,2));
% y_gauss = zeros(nQuad,size(IEN,2));
% 
% % elastic strain of gauss point
% strainvector_gauss_x = zeros(nQuad,size(IEN,2));
% strainvector_gauss_y = zeros(nQuad,size(IEN,2));
% strainvector_gauss_xy = zeros(nQuad,size(IEN,2));
% 
% % eigen strain
% global eigen_strain_matrix
% eigen_strain_matrix= zeros(nQuad,size(IEN,2));

% global coord_psi_matrix
% coord_psi_matrix= zeros(nQuad,size(IEN,2));

% stress of gauss point
% global strain_x strain_y strain_z strain_xy   
% strain_x= zeros(nQuad,size(IEN,2));
% strain_y= zeros(nQuad,size(IEN,2));
% strain_z= zeros(nQuad,size(IEN,2));
% strain_xy= zeros(nQuad,size(IEN,2));
% 
% global lambda_max global_lambda_max
% lambda_max = zeros(nQuad,size(IEN,2));
% global_lambda_max = zeros(nQuad,size(IEN,2));

% global pl_energy_all pl_energy_all_old
% pl_energy_all = zeros(nQuad,size(IEN,2));
% pl_energy_all_old= zeros(nQuad,size(IEN,2));

% global H_gauss_record  
% %alpha_bar_gauss_record
% H_gauss_record = zeros(nQuad * size(IEN,2),1);  
% % alpha_bar_gauss_record = zeros(nQuad * size(IEN,2),1);

global temperature 
temperature=zeros(nQuad,size(IEN,2));

global Sol_T_record PI_ud
Sol_T_record = zeros(length(Sol_T),1);
PI_ud = 0;


% RuR = zeros(nNodes,1);
% RuCoordR = zeros(nNodes,1);
% RuCoord_phi = zeros(nNodes,1);
% Ru_record = zeros(maxloadstep,1);

% H_gauss = zeros(nQuad * size(IEN,2),500);%500 gives randomly
% alpha_bar_gauss = zeros(nQuad * size(IEN,2),500);

% nabla_d_last = [zeros(1,size(Coord,2));ones(1,size(Coord,2))];

% plastic
% dev_ep_LastLoadStep = zeros(6*nQuad,size(IEN,2)); %nQuad by nElements, eff plastic strain
stressvec = zeros(3*nQuad,size(IEN,2)); %nQuad by nElements, stress vector
% plstrain_old_tensor = zeros(6*nQuad,size(IEN,2));%plastic strain
% IValpha_LastLoadStep = zeros(6 * nQuad,size(IEN,2)); %nQuad on each nElements,the size of IValpha is 6 on each integration point
% plStrainvector_LastLoadStep = zeros(nNodesElement*nDoF,nElements);      %plastic strain


psi_plus_vector = zeros(maxloadstep,2);
% pure_psi_plus_vector = zeros(maxloadstep,2);
% c_length_vector = zeros(maxloadstep,2);
% d_matrix = zeros(nQuad,size(IEN,2),maxloadstep);
loadpre_vector = zeros(maxloadstep,1);
% Ru_average_vector = zeros(maxloadstep,2);
% Reaction_average_vector = zeros(maxloadstep,2);
c_length_recorder = zeros(500,5);   % record  the length of crack

%variables for force control
force = zeros(maxloadstep,1);
u_node = zeros(maxloadstep,1);
c_length = zeros(maxloadstep,1);%to represent the crack length


sum_psi_plus_rec = zeros(maxloadstep,1);


load_pre =0;   % Current load time step
sum_crack_RR = 0;
% load_max = 1; %close cycle
% load_min = -0.4;

% pic_num = 1;
deltau =1e-4; %load time step %initial loadstep should be smaller than load_max
load_control = deltau * 5;  %control load step
step_no = 0;
Nt_iter_no1 = 1;
Nt_iter_no2 = 1;
Nt_iter_noT = 1;    %temperature
max_iter = 1000;
sumRu = zeros(maxloadstep,2);
%point = 1;
% zero_one_vector = zeros(5000,1);
% loadrecord = zeros(10000,2);
% for iter = 1 : length(zero_one_vector)
%     if mod(iter,2)
%         zero_one_vector(iter) = 1;
%     end
% end
conver_vector = zeros(max_iter,1);
converge=1;
% Loadstep_recorder = zeros(5,2);
% falphabar_recorder = zeros(5,1000);
% r1 = 0;

% R_tol_u = 1e-6;
Sol_u_LastLoadStep = Sol_u;
Sol_d_LastLoadStep = Sol_d;
Sol_T_record = Sol_T;
flag1 = 's求解';
flag2 = 's求解';
m = 0.95;
% n = 1e-2;
rho_ref = 0.1;
Ru_tol = 1e-6;Rd_tol = 1e-6;RT_tol = 1e-6; % 配合R_ref的 （放在这里是因为后面可能会改动所以，每一步都要刷新这个）
Ru_ano = 1e-4;Rd_ano = 1e-15;RT_ano = 1e-3; % 绝对值的容差
%%
rec_RTs = zeros(1,3);
rec_RTsc = zeros(1,3);
rec_Rus = zeros(1,3);
rec_Rusc = zeros(1,3);
 % load initial0ud.mat
record = zeros(maxloadstep,9);
% Ru_ref = 6.26e4;Rd_ref = 0.005;RT_ref = 7.6e4;
% Ru_tol = 1e-3;Rd_tol = 1e-3;RT_tol = 1e-6;
%         load Copy_of_1.11.mat;%contiue the computation

% max_iter = 100;
% period=1;%input 1/8 time of one cycle loading
%for no = 1:maxloadstep
for no = 1:96
    iterud = 0; % 为了记录当前加载步下ud之间的迭代次数
    open_eta1 = 0;
    open_eta2 = 0;
    runtime=0;
    tic
%     if(point == 480)
%         pause;
%     end   
    step_no = step_no + 1;  
    fprintf('%g th load steps\n',step_no);
    
    converged = false;

    load_old = load_pre;
%    sum_crack_RR_old = sum_crack_RR;
%     if load_old>=1.5e-3
%     deltau = 5e-6;
%     end
%     %
%     fixed = 1;
%     semi_period_time = 1;
%     quarter_period_time = 0.5;
%     quarter_cycle = 0.5*fix(load_pre/quarter_period_time) + 0.5;
%     semi_cycle = fix(load_pre/semi_period_time) + 1;
%     dt_max = semi_cycle - load_pre;
%     dt_max2 = quarter_cycle - load_pre;
%     if dt_max2 >= deltau
        load_pre = load_pre + deltau;

%       temperature=4*load_pre+20;
%     else
%         load_pre = load_pre + dt_max2;
%     end
%     %
    loadrecord(step_no,:) = [step_no load_pre];
    fprintf('load_pr esent is %f ,deltau is %f\n',load_pre,deltau);

   
%____________________temperature____________________________
    Sol_T(BCIndices_T) = BCVal_T(BCIndices_T);   
    Sol_T_record = Sol_T;
    switch iterationType
        case 'monolithic'
            Sol(BCIndices) = BCVal(BCIndices) * load_pre;
            %load steps are to ease the solution of nonlinear problems.
            % It can be taken as 1 or [0,1], or something like 0:0.1:1
            % The body force and traction are scaled by the same factor
        case 'staggering'
            Sol_u(BCIndices_u) = BCVal_u(BCIndices_u);
            % Sol_d(BCIndices_d) = BCVal_d(BCIndices_d) * load_steps(step_no);
            Sol_d(BCIndices_d) = BCVal_d(BCIndices_d);
    end
    scu_times = 0;
    ntu_times = 0;
    scd_times = 0;
    ntd_times = 0;
    scT_times = 0;
    ntT_times = 0;
    exa_scutimes = 0;
    
    
    
    for iter_noall=1:max_iter
        %for i = 1:2
        
        nt_T_converged = false;
        Ru = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
        Rd = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
        RT = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
            constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
        Ru_active = Ru(active_indices_u);
        Rd_active = Rd(active_indices_d);
        RT_active = RT(active_indices_T);

        if strcmp(flag1,'s求解')
            rec_RTs = [rec_RTs,norm(RT_active)];rec_RTs = rec_RTs(:,2:4);
            if ((rec_RTs(3)/rec_RTs(2)>rho_ref)||(rec_RTs(2)/rec_RTs(1)>rho_ref))&&(rec_RTs(1)*rec_RTs(2)*rec_RTs(3)~=0)
                flag1='sc法';
            end
        else
            rec_RTsc = [rec_RTsc,norm(RT_active)];rec_RTsc = rec_RTsc(:,2:4);
            if ((rec_RTsc(3)/rec_RTsc(2)>1)||(rec_RTsc(2)/rec_RTsc(1)>1))&&(rec_RTsc(1)*rec_RTsc(2)*rec_RTsc(3)~=0)
                open_eta1 = 1;
            end
        end
        fprintf('L2-norm of Ru = %g, L2-norm of Rd = %g, L2-norm of RT = %g\n',norm(Ru_active),norm(Rd_active),norm(RT_active));
        if iter_noall==1
            RT_ref = norm(RT_active);
            for nt_T=1:max_iter
                [RT,~] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
                RT_active = RT(active_indices_T);
                fprintf('L2-norm of RT = %g\n',norm(RT_active));
                %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
                if (norm(RT_active)<= RT_tol * RT_ref)
                    break;
                end
                [~,KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
                KT_active = KT(active_indices_T, active_indices_T);
                Sol_T(active_indices_T) = Sol_T(active_indices_T) - KT_active\RT_active;
                %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
            end
            Ru = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
            Rd = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
            Ru_active = Ru(active_indices_u);
            Rd_active = Rd(active_indices_d);
            Ru_ref = norm(Ru_active);
            Rd_ref = norm(Rd_active);
        end
        if (norm(Ru_active)<= Ru_tol * Ru_ref) &&  (norm(Rd_active)<= Rd_tol * Rd_ref) && (norm(RT_active)<= RT_tol * RT_ref)
            converged = true;
            break;
        end
        if iter_noall ~= 1
            if strcmp(flag1,'s求解')
                for nt_T=1:max_iter
                    [RT,~] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
                    RT_active = RT(active_indices_T);
                    fprintf('L2-norm of RT = %g\n',norm(RT_active));
                    %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
                    if (norm(RT_active)<= RT_tol * RT_ref)
                        break;
                    end
                    [~,KT] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'T');
                    KT_active = KT(active_indices_T, active_indices_T);
                    Sol_T(active_indices_T) = Sol_T(active_indices_T) - KT_active\RT_active;
                    %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
                end
            elseif strcmp(flag1,'sc法')
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
                eta1 = 1;% 默认等于1
                if open_eta1 == 1
                    rho1 = abs(eigs([sparse(size(KT_active, 1), size(KT_active, 2)), [KTu_active, KTd_active]; ...
                        sparse(size([KuT_active; KdT_active], 1), size([KuT_active; KdT_active], 2)), sparse(size(Ku_active, 1) + size(Kd_active, 1), size(Ku_active, 2) + size(Kd_active, 2))], ...
                        [KT_active, sparse(size([KTu_active, KTd_active],1),size([KTu_active, KTd_active],2)); ...
                        [KuT_active; KdT_active], [Ku_active, Kud_active; Kdu_active, Kd_active]], 1));
                    if rho1<1
                        eta1 = 1;
                        if rho1<rho_ref
                            % flag1 = 's求解';
                        end
                    else
                        eta1 = m/sqrt(rho1);
                        fprintf('eta1 = %d\n',eta1)
                    end
                end
                K_sc = KT_active - eta1^2*[KTu_active,KTd_active]/[Ku_active,Kud_active; Kdu_active,Kd_active]*[KuT_active;KdT_active];
                R_sc = RT_active - eta1*[KTu_active,KTd_active]/[Ku_active,Kud_active; Kdu_active,Kd_active]*[Ru_active;Rd_active];
                Sol_T(active_indices_T) = Sol_T(active_indices_T) - K_sc\R_sc;
            end
        end
        for iter_ud = 1:max_iter
            % 接下来更新ud
            Rd = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
            Ru = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
            Rd_active = Rd(active_indices_d);
            Ru_active = Ru(active_indices_u);
            fprintf('L2-norm of Ru = %g, L2-norm of Rd = %g\n',norm(Ru_active),norm(Rd_active));
            
            if (norm(Ru_active)<= Ru_tol * Ru_ref) && (norm(Rd_active)<= Rd_tol * Rd_ref) 
                    converged = true;
                    break;
            end
            iterud = iterud + 1; %记录在此加载步下ud之间迭代了多少次
            if strcmp(flag2,'s求解')
                rec_Rus = [rec_Rus,norm(Ru_active)];rec_Rus = rec_Rus(:,2:4);
                if (rec_Rus(3)/rec_Rus(2)>rho_ref)&&(rec_Rus(2)/rec_Rus(1)>rho_ref)&&(rec_Rus(1)*rec_Rus(2)*rec_Rus(3)~=0)
                    flag2='sc法';
                end
            else
                rec_Rusc = [rec_Rusc,norm(Ru_active)];rec_Rusc = rec_Rusc(:,2:4);
                if ((rec_Rusc(3)/rec_Rusc(2)>1)||(rec_Rusc(2)/rec_Rusc(1)>1))&&(rec_Rusc(1)*rec_Rusc(2)*rec_Rusc(3)~=0)
                    open_eta2 = 1;
                end
            end
            if strcmp(flag2,'s求解')% d用什么
                for nt_u=1:max_iter
                    [Ru,~] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
                    Ru_active = Ru(active_indices_u);
                    fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
                    %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
                    if (norm(Ru_active)<= Ru_tol * Ru_ref)
                        break;
                    end
                    [~,Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                        constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
                    Ku_active = Ku(active_indices_u, active_indices_u);
                    Sol_u(active_indices_u) = Sol_u(active_indices_u) - Ku_active\Ru_active;
                    %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
                end
            elseif strcmp(flag2,'sc法')
                [~,Kud] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'ud');
                [~,Kdu] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'du');
                [~,Ku] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'u');
                [~,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
                Kd_active = Kd(active_indices_d, active_indices_d);
                Kdu_active = Kdu(active_indices_d, active_indices_u);
                Kud_active = Kud(active_indices_u, active_indices_d);
                Ku_active = Ku(active_indices_u, active_indices_u);
                eta2 = 1;
                if open_eta2 == 1
                    rho2 = abs(eigs([sparse(size(Ku_active, 1), size(Ku_active, 2)), Kud_active; ...
                        sparse(size(Kdu_active, 1), size(Kdu_active, 2)), sparse(size(Kd_active, 1), size(Kd_active, 2))], ...
                        [Ku_active, sparse(size(Kud_active, 1), size(Kud_active, 2)); ...
                        Kdu_active, Kd_active], 1));
                   if rho2<1
                       eta2 = 1;
                       if rho2<rho_ref
                           flag2 = 's求解';
                       end
                   else
                       eta2 = m/sqrt(rho2);
                       fprintf('eta2 = %d\n',eta2)
                   end
                end
                K_sc = Ku_active - eta2^2*Kud_active/Kd_active*Kdu_active; % dT的舒尔补
                R_sc = Ru_active - eta2*Kud_active/Kd_active*Rd_active;
                Sol_u(active_indices_u) = Sol_u(active_indices_u) - K_sc\R_sc;

            end

            for nt_d=1:max_iter
                [Rd,~] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
                Rd_active = Rd(active_indices_d);
                fprintf('L2-norm of Rd = %g\n',norm(Rd_active));
                %     if norm(RT_active)<= R_tol * sqrt(length(RT_active))/102.6
                if (norm(Rd_active)<Rd_tol * Rd_ref)
                    break;
                end
                [~,Kd] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
                    constitutive,body_force,traction, load_pre,step_no,nQuad,gc_vec,'d');
                Kd_active = Kd(active_indices_d, active_indices_d);
                Sol_d(active_indices_d) = Sol_d(active_indices_d) - Kd_active\Rd_active;
                %d_times = d_times + 1;%Increment the number of displacement field iterations by one.
            end

        end    

    end


    if iter_noall==max_iter
        converge=0;
        break;
    end

    if converge==0
        fprintf('There are too many RT Newton iterations!\n');
        break;
    end

    %             X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
    %             Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
    %             Sol_T_c = Sol_T - 273;
    %             answer = Sol_T_c';
    %             T_elem = [answer(IEN(1,:)); answer(IEN(2,:));answer(IEN(3,:))];
    %
    %             figure(1)
    %             colormap('jet');
%             patch(X,Y,T_elem);
%             axis equal;
%             axis off;
%              shading interp;
%             colorbar;
%             caxis([20 300]);
%             title(sprintf('temperature under loading step %d ', step_no));
%             h1=gcf;
%             saveas(h1,['C:\Users\吕冰\Desktop\温度场50步\',num2str(step_no),'.jpg']);
% end

% Sol_T_LastLoadStep = Sol_T;

 %__________________________________________________________

     if (converged)
        converge_laststep = true;
        runtime = toc
        %update the history variables
        Psi_plus_old = Psi_plus_rec;
        record(no,1) = load_pre;
        record(no,2) = scu_times;
        record(no,3) = ntu_times;
        record(no,4) = scd_times;
        record(no,5) = ntd_times;
        record(no,6) = scT_times;
        record(no,7) = iterud;
        record(no,8) = iter_noall;
        record(no,9) = runtime;
        Sol_u_LastLoadStep = Sol_u;
        Sol_d_LastLoadStep = Sol_d;
        Sol_T_LastLoadStep = Sol_T;

        %% if converged, output the results
%         load test.mat
        F = Ru(BCIndices_u);
        sum(abs(F));
        if strcmp(iterationType,'staggering')
            Sol(u_indices) = Sol_u;
            Sol(d_indices) = Sol_d;
            Sol(T_indices) = Sol_T;
        end
        
        switch problemType
            case 'temperature'
                Answer=reshape(Sol, 4, []);
                Displace=Answer(1:2,:);
%                 Displace_x=Answer(1,:);
%                 Displace_y=Answer(2,:);
                D=Answer(3,:);
            case 'PhaseField'
                Answer=reshape(Sol, 3, []);
                Displace=Answer(1:2,:);
                D=Answer(3,:);
            case 'Poisson'
                D=Sol';
        end
     
        
        X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
        Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
        D_elem = [D(IEN(1,:)); D(IEN(2,:)); D(IEN(3,:))];
        
%         
        figure(2)
        colormap('jet');
        XYZ = patch(X,Y,D_elem);
        axis equal;
        axis off;
        shading interp;
        colorbar;
        caxis([0 1]);
        title(sprintf('The phase field under loading step %d ', step_no));
%         h1=gcf;
% %         if rem(step_no,5)==0
%               saveas(h1,['C:\Users\吕冰\Desktop\519\','nt',num2str(step_no),'.jpg']);
%         end
% 
%         figure(3)
%         Coord_deformation=Coord+Displace;
%         X = [Coord_deformation(1,IEN(1,:)); Coord_deformation(1,IEN(2,:)); Coord_deformation(1,IEN(3,:))];
%         Y = [Coord_deformation(2,IEN(1,:)); Coord_deformation(2,IEN(2,:)); Coord_deformation(2,IEN(3,:))];
%         ux = [Displace(1, IEN(1,:)); Displace(1, IEN(2,:)); Displace(1, IEN(3,:))];
%         uy = [Displace(2, IEN(1,:)); Displace(2, IEN(2,:)); Displace(2, IEN(3,:))];
%         patch(X,Y,ux);
%         axis equal;
%         axis off;
%         shading interp;
%         colorbar;
%         title(sprintf('The deformation under loading step %d ', step_no));
% 
%         figure(4)
%         X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
%         Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
%         ux = [Displace(1, IEN(1,:)); Displace(1, IEN(2,:)); Displace(1, IEN(3,:))];
%         uy = [Displace(2, IEN(1,:)); Displace(2, IEN(2,:)); Displace(2, IEN(3,:))];
%         patch(X,Y,ux)
%         axis equal;
%         axis off;
%         shading interp;
%         colorbar;
%         title(sprintf('The x-displacement under loading step %d ', step_no));
%         h2=gcf;
%         saveas(h2,['C:\Users\吕冰\Desktop\519\','x位移',num2str(step_no),'.jpg']);
% 
% 
%         figure(5)
%         patch(X,Y,uy)
%         axis equal;
%         axis off;
%         shading interp;
%         colorbar;
%         title(sprintf('The y-displacement under loading step %d ', step_no));
%         h3=gcf;
%         saveas(h3,['C:\Users\吕冰\Desktop\4e-6\y位移\',num2str(step_no),'png']);
% 
%          figure(6)
%         s11 = [stressvec_old(1,:); stressvec_old(4,:); stressvec_old(7,:)];
%         s22 = [stressvec_old(2,:); stressvec_old(5,:); stressvec_old(8,:)];
%         s12 = [stressvec_old(3,:); stressvec_old(6,:); stressvec_old(9,:)];
%         patch(X,Y,s22);     
%         shading interp;       
%         colorbar;            
%         axis equal;          
%         axis equal;
%         axis off;
%         title(sprintf('/sigma_{y} under loading step %d ', step_no)); 
%         h4=gcf;
%         saveas(h4,['C:\Users\吕冰\Desktop\4e-6\y应力\',num2str(step_no),'png']);
% 
%         figure(7)
%         patch(X,Y,global_lambda_max)
%         axis equal;
%         axis off;
%         shading interp;
%         colorbar;
%         title(['\sigma_{max} under loading step', num2str(step_no)]);
%         h3 = gcf;
%         saveas(h3,['C:\Users\吕冰\Desktop\519\','lambda',num2str(step_no),'.jpg']);
        

%         figure(8)
%         patch(X,Y,lambda_max)
%         axis equal;
%         axis off;
%         shading interp;
%         colorbar;
%         title(sprintf('\sigma_{max} under loading step %d ', step_no));
%         h5=gcf;
%         saveas(h5,['C:\Users\吕冰\Desktop\lambda\',num2str(step_no),'png']);

%         figure(8)
%         patch(X,Y,s11)
%         axis equal;
%         axis off;
%         shading interp;
%         colorbar;
%         title(sprintf('sigma_{x} under loading step %d ', step_no));
%         h6=gcf;
%         saveas(h6,['C:\Users\吕冰\Desktop\4e-6\x应力\',num2str(step_no),'png']);

            fprintf('maximum d = %f\n',max(max(Sol_d)))

%         principal_stress_max(step_no) = max(max(lambda_max));
%         time(step_no) = load_pre*10;
%         figure(9)
%         plot(time,principal_stress_max);
%         
%         figure(10)
%         plot(time,nc_PI);
%         ylabel('energy[J]');
%         xlabel('time[ms]');
%         grid on

%drawing dynamic picture
%         F=getframe(gcf);
%         I=frame2im(F);
%         [I,map]=rgb2ind(I,256);
% 
%         if pic_num0 == 1
%             imwrite(I,map,'d.gif','gif', 'Loopcount',inf,'DelayTime',0.1);
%         else
%             imwrite(I,map,'d.gif','gif','WriteMode','append','DelayTime',0.1);
%         end
%         pic_num0 = pic_num0 + 1;

%         title(sprintf('The d under loading step %d ', step_no));
      

%         period_time=2;
%         cycle=fix(load_pre/period_time)+1;
       

%         ep_elem = [ep_new_vec(1,:); ep_new_vec(2,:); ep_new_vec(3,:)];    %%plStrainvector in yy direction

        
        
%          c_length(step_no) = sum(crack_RR);
%         if mod(load_pre,period) == 0
%             c_length_recorder(fix(load_pre/period)+1,:) = [fix(load_pre/period)+1 sum(crack_RR) sum(sum(Psi_plus_old)) sum(sum(Psi_plus_old_d)) sum(sum(driving_force))];
%         end
%         
% 
%          fprintf('c_length(%g) = %g under t= %g\n',step_no,c_length(step_no),load_pre);


  


%         sum_alpha(step_no) = sum(sum(alpha_bar_new));
%         fprintf('sum_alpha(%g) = %g under t= %g\n',step_no,sum_alpha(step_no),load_pre);
%         norm_alpha(step_no) = norm(alpha_bar_new);
%         fprintf('norm_alpha(%g) = %g under t= %g\n',step_no,norm_alpha(step_no),load_pre);
%         sum_psi_plus_rec(step_no) = sum(sum(Psi_plus_rec));
%         fprintf('sum_psi_plus_rec(%g) = %g under t= %g\n',step_no,sum_psi_plus_rec(step_no),load_pre);

        loadpre_vector(step_no,1) = load_pre;
        %alpha_matrix(:,:,step_no) = alpha_bar_new;
        d_matrix(:,:,step_no) = D_elem;


%         figure(9)
%         norm_alpha_vector(step_no,2) = load_pre;
%         norm_alpha_vector(step_no,1) = norm_alpha(step_no);
%         scatter(norm_alpha_vector(:,2),norm_alpha_vector(:,1),'filled');
%         title(sprintf('The norm of alpha under t= %.6f ', load_pre));
%         xlabel('t') ;
%         ylabel('norm of alpha') ;
%         grid on;



         psi_plus_vector(step_no,1) = sum(sum(Psi_plus_old));
         psi_plus_vector(step_no,2) = load_pre;

        
%         subplot(2,2,3);
%         alpha_bar_vector(step_no,1) =  sum(sum(alpha_bar_old));
%         alpha_bar_vector(step_no,2) = load_pre;
%         scatter(alpha_bar_vector(:,2),alpha_bar_vector(:,1));
%         title(sprintf('The alpha_bar_sum under loading step %d ', step_no));
%         xlabel('defomation') ;
%         ylabel('alpha_bar_sum') ;
        



%         pic_num=pic_num+1;
        

%         if load_pre==0.5  
%         filename = 'test_quater.mat';
%         save(filename);
%         break;
%         end
%         if load_pre <= 400.5  
%         filename = 'test.mat';
%         save(filename);
%         end
%         if load_pre==400.5
%         break;
%         end

%         if load_pre<=100000
% %             filename = '1.11.mat';
% %             filename = 'test5.mat';%第12步 温度和位移交换位置
% %              filename = 'test.mat';%温度耦合相场 2e-5
% %              filename = 'test6.mat';
%             save(filename);
%         end

%         if mod(fix(load_pre/period)+1,5001) == 0
%         filename50x = ['test',num2str(fix(load_pre/period)+1,'.mat')];
%         save(filename50x);
%         end

    else
        load_pre = load_old;
        disp('//////////////////////////////////////////');
        fprintf('load_pre returned to %f\n',load_pre)
        deltau = deltau/2;
        if deltau > abs(load_control)/5
            deltau = abs(load_control)/5;
        end
        Sol_u = Sol_u_LastLoadStep;
        Sol_d = Sol_d_LastLoadStep;
        disp('reduce the loadstep, and begin from the last step');
        disp('//////////////////////////////////////////');
     end
     saveas(XYZ,['Thermo8467ele/Phasefield-' strrep(num2str(load_pre, '%.5f'),'.','_') '.jpg']);
     clear K_sc R_sc Ku Kd Kud Kdu KuT KTu KdT KTd
     save(['Thermo8467ele/PF-' strrep(num2str(load_pre, '%.5f'),'.','_') '.mat']);
end

% save ('Ru_average.mat', 'Ru_record');
%---------4/17-------------------
% save in 4.17, Using an unmodified algorithm 
% exp(-3)
% plot the energy vs time
% 9th
%----------------------------------

%---------5-4----------------------
% R_tol = 1e-7
% draw the energy picture
% ----------------------------------
%舒尔补先scud，scdu如果不行就牛顿法求至收敛，然后交替法求T
%
