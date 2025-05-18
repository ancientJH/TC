%%
clc
clear
elementType = 'P12D';
ElementType = 'P12D';% 'P12D' means 2D triangular element,'P22D' means 2D triangular element with 2-order shape function
problemType = 'PhaseField';         % Other options: 'Elasticity''PhaseField'

[nDim, nDoF, nNodes, nElements, nNodesElement, Coord, nEquations, ...
    DirichletBCs,IEN, LM, constitutive, body_force, traction] = ...
    ProblemDefinition(elementType, problemType);

%% mesh
figure(1)
X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
%D_elem = ones(size(X));
D_elem = linspace(1,size(IEN,2),size(IEN,2));

% Patch

hold on
XYZ=patch(X,Y,D_elem);
axis equal;
axis off;
colorbar;


%generate gauss points on each element
% gauss_on_ele_y = zeros(size(IEN));
% gauss_on_ele_x = zeros(size(IEN));
switch ElementType
    case 'P12D'
        nQuad = 3; nQuadBdry = 2;
    case 'Q12D'
        nQuad = 4; nQuadBdry = 2;
end

%--------------------------------------------------------------------------
% Phase 2: Solve with load stepping
%--------------------------------------------------------------------------
%load_steps = -load_steps;
% figure(2)
% plot(load_steps);
% hold on;
% plot(cyc_record);
% title(sprintf('The load steps'));
% xlabel('time step');
% ylabel('load value');



maxloadstep = 100000;
Sol = zeros(nEquations, 1);
[bc_check,BCIndices, BCVal,gc_vec] = DirichletBCs(Coord);

% subplot(2,1,2);%下面可以注释掉
% X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
% Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
% gc_node = [gc_vec(IEN(1,:));gc_vec(IEN(2,:));gc_vec(IEN(3,:))];
% bc_node_d = [bc_check(3,IEN(1,:)); bc_check(3,IEN(2,:)); bc_check(3,IEN(3,:))];
% hold on
% XYgc=patch(X,Y,bc_node_d);
% axis equal;
% axis off;
% colorbar;

if strcmp(problemType, 'Poisson') && sum(BCIndices) == 0 % fix the value of any node in case of no Dirichlet BC
    BCIndices(1) = true;
    BCVal(1) = 0;
end
active_indices = ~BCIndices;   %Not on the Dirichlet boundary, i.e., free boundary
R_tol = 1e-14;    %----------->absolute tolerance of the residual
% R_tol = 1e-40; 
iterationType = 'staggering'; % Possible values: 'monolithic' 整体法or 'staggering'交替法

if strcmp(iterationType, 'staggering')
    if strcmp(problemType, 'Elasticity') || strcmp(problemType, 'Poisson')
        iterationType = 'monolithic';
    else
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
ll=0;
% pic_num0 = 1;
% pic_num1 = 1;
% pic_num2 = 1;
% pic_num3 = 1;
global Psi_plus_rec Psi_plus_old Psi_plus_old_d Psi_plus_rec_d Pure_Psi_plus driving_force
Psi_plus_old = zeros(nQuad,size(IEN,2)); %4 is the number of gauss points.size(IEN,2) is the size of the elements
Psi_plus_rec = zeros(nQuad,size(IEN,2));
Psi_plus_old_d = zeros(nQuad,size(IEN,2));
Psi_plus_rec_d = zeros(nQuad,size(IEN,2));
Pure_Psi_plus = zeros(nQuad,size(IEN,2));
driving_force = zeros(nQuad,size(IEN,2));

%add function f(alpha_bar t) in the artical"A framework to model the fatigue behavior of brittle materials based
%on a variational phase-field approach"
global alpha_bar_old alpha_bar_new f_alpha_bar
alpha_bar_old = zeros(nQuad,size(IEN,2));%4 is the number of gauss points.
alpha_bar_new = zeros(nQuad,size(IEN,2));
f_alpha_bar = zeros(nQuad,size(IEN,2));

global stressvector_gauss stressvector_gauss_r stressvector_gauss_z Hydrostatic_stress strainvector_theta_gauss 
stressvector_gauss = zeros(nQuad,size(IEN,2));%4 is the number of gauss points.It is the y direction of the gauss points, you can change the direction by
stressvector_gauss_r = zeros(nQuad,size(IEN,2));
stressvector_gauss_z = zeros(nQuad,size(IEN,2));
Hydrostatic_stress = zeros(nQuad,size(IEN,2));
strainvector_theta_gauss = zeros(nQuad,size(IEN,2));
global R_gauss rr_gauss zz_gauss strainvector_r_gauss strainvector_z_gauss 
R_gauss = zeros(nQuad,size(IEN,2));
rr_gauss = zeros(nQuad,size(IEN,2));
zz_gauss = zeros(nQuad,size(IEN,2));
strainvector_r_gauss = zeros(nQuad,size(IEN,2));
strainvector_z_gauss = zeros(nQuad,size(IEN,2));
%changing the stressvector_gauss(iQuad,ielem) = StressVector(2) in PhaseFieldElement.m
% 
global coord_stressvector_guass_R coord_stressvector_guass_theta coord_stressvector_guass_psi
coord_stressvector_guass_R = zeros(nQuad,size(IEN,2));
coord_stressvector_guass_theta = zeros(nQuad,size(IEN,2));
coord_stressvector_guass_psi = zeros(nQuad,size(IEN,2));
% 
global eigen_strain_matrix
eigen_strain_matrix= zeros(nQuad,size(IEN,2));

global strain_r strain_theta strain_z strain_rz
strain_r= zeros(nQuad,size(IEN,2));
strain_theta= zeros(nQuad,size(IEN,2));
strain_z= zeros(nQuad,size(IEN,2));
strain_rz= zeros(nQuad,size(IEN,2));


global coord_psi_matrix
coord_psi_matrix= zeros(nQuad,size(IEN,2));
% 
global lambda_max lambda_1 lambda_2 lambda_3  lambda_max_guass lambda_max_node 
lambda_max = zeros(nQuad,size(IEN,2));
lambda_1 = zeros(nQuad,size(IEN,2));
lambda_2 = zeros(nQuad,size(IEN,2));
lambda_3 = zeros(nQuad,size(IEN,2));
lambda_max_guass = zeros(size(IEN,2));
lambda_max_node = 0;



global pl_energy_all pl_energy_all_old
pl_energy_all = zeros(nQuad,size(IEN,2));
pl_energy_all_old= zeros(nQuad,size(IEN,2));
global H_gauss_record alpha_bar_gauss_record
H_gauss_record = zeros(nQuad * size(IEN,2),1);
alpha_bar_gauss_record = zeros(nQuad * size(IEN,2),1);

RuR = zeros(nNodes,1);
RuCoordR = zeros(nNodes,1);
RuCoord_phi = zeros(nNodes,1);
Ru_record = zeros(maxloadstep,1);

H_gauss = zeros(nQuad * size(IEN,2),500);%500 gives randomly
alpha_bar_gauss = zeros(nQuad * size(IEN,2),500);

nabla_d_last = [zeros(1,size(Coord,2));ones(1,size(Coord,2))];

dev_ep_LastLoadStep = zeros(6*nQuad,size(IEN,2)); %nQuad by nElements, eff plastic strain
stressvec = zeros(3*nQuad,size(IEN,2)); %nQuad by nElements, stress vector
plstrain_old_tensor = zeros(6*nQuad,size(IEN,2));%plastic strain
IValpha_LastLoadStep = zeros(6 * nQuad,size(IEN,2)); %nQuad on each nElements,the size of IValpha is 6 on each integration point
plStrainvector_LastLoadStep = zeros(nNodesElement*nDoF,nElements);      %plastic strain

ep_vector = zeros(size(IEN));
alpha_bar_vector = zeros(maxloadstep,2);
f_alpha_bar_vector = zeros(maxloadstep,2);
psi_plus_vector = zeros(maxloadstep,2);
pure_psi_plus_vector = zeros(maxloadstep,2);
c_length_vector = zeros(maxloadstep,2);
sum_alpha_vector = zeros(maxloadstep,2);%
norm_alpha_vector = zeros(maxloadstep,2);%疲劳
%alpha_matrix = zeros(nQuad,size(IEN,2),maxloadstep);
%d_matrix = zeros(nQuad,size(IEN,2),maxloadstep);
loadpre_vector = zeros(maxloadstep,1);
Ru_average_vector = zeros(maxloadstep,2);
Reaction_average_vector = zeros(maxloadstep,2);
c_length_recorder = zeros(500,5);%记录裂纹长度
%variables for force control
force = zeros(maxloadstep,1);
u_node = zeros(maxloadstep,1);
c_length = zeros(maxloadstep,1);%to represent the crack length
sum_alpha = zeros(maxloadstep,1);
sum_alpha = zeros(maxloadstep,1);
norm_alpha = zeros(maxloadstep,1);
sum_psi_plus_rec = zeros(maxloadstep,1);
eigenstrain_max=0.0225504;
% eigenstrain_max=0.01;
radius=0.0005;
load_pre = -0.25;
sum_crack_RR = 0;
% load_max = 1; %close cycle
% load_min = -0.4;
pic_num = 1;
deltau = 0.25;%initial loadstep should be smaller than load_max时间
load_control = deltau * 5;
iter_no = 1;
step_no = 0;
Nt_iter_no1 = 1;
Nt_iter_no2 = 1;
max_iter = 100;
max_iter1 = 100;
max_iter2 = 100;
sumRu = zeros(maxloadstep,2);
%point = 1;
zero_one_vector = zeros(5000,1);
loadrecord = zeros(10000,2);
Ruactive_vector = zeros(max_iter1,1);
for iter = 1 : length(zero_one_vector)
    if mod(iter,2)
        zero_one_vector(iter) = 1;
    end
end
conver_vector = zeros(max_iter,1);
Loadstep_recorder = zeros(5,2);
falphabar_recorder = zeros(5,1000);
r1 = 0;
% % R_tol = 1e-6;
R_tol = 1e-6;
Sol_u_LastLoadStep = Sol_u;
Sol_d_LastLoadStep = Sol_d;
%%
load test.mat;%contiue the computation
% load test_quater.mat;
%运行节！！！
period=1;%input 1/8 time of one cycle loading
for no = 1 : maxloadstep
%     if(point == 480)
%         pause;
%     end   
    step_no = step_no + 1;
    fprintf('%g th load steps\n',step_no);
    
    converged = false;

    load_old = load_pre;
    sum_crack_RR_old = sum_crack_RR;
    
    fixed = 1;
    semi_period_time = 1;
    quarter_period_time = 0.5;
    quarter_cycle = 0.5*fix(load_pre/quarter_period_time) + 0.5;
    semi_cycle = fix(load_pre/semi_period_time) + 1;
    dt_max = semi_cycle - load_pre;
    dt_max2 = quarter_cycle - load_pre;
    if dt_max2 >= deltau
        load_pre = load_pre + deltau;
    else
        load_pre = load_pre + dt_max2;
    end

    loadrecord(step_no,:) = [step_no load_pre];
    fprintf('load_present is %f ,deltau is %f\n',load_pre,deltau);
    
    switch iterationType
        case 'monolithic'
            Sol(BCIndices) = BCVal(BCIndices) * load_pre;
            %load steps are to ease the solution of nonlinear problems.
            % It can be taken as 1 or [0,1], or something like 0:0.1:1
            % The body force and traction are scaled by the same factor
        case 'staggering'
            Sol_u(BCIndices_u) = BCVal_u(BCIndices_u) * load_pre;
            % Sol_d(BCIndices_d) = BCVal_d(BCIndices_d) * load_steps(step_no);
            Sol_d(BCIndices_d) = BCVal_d(BCIndices_d);
    end
    
    for iter_no = 1:max_iter  %-------->Stepping for Newton iteration
        fprintf('%g th iter_no\n',iter_no);
        switch iterationType
            case 'monolithic'
                R = COHESIVEAssemble(Sol, Coord, IEN, LM, elementType, problemType, ...   %R is the global residual vector
                    constitutive, body_force, traction, load_pre);
                R_active = R(active_indices);
                fprintf('L2-norm of residual = %g\n',norm(R_active));
                if norm(R_active) <= R_tol * sqrt(length(R_active))
                    converged = true;
                    break;
                end
                
                [~, K] = COHESIVEAssemble(Sol, Coord, IEN, LM, elementType, problemType, ...
                    constitutive, body_force, traction, load_pre);
                K_active = K(active_indices, active_indices);
                Sol(active_indices) = Sol(active_indices) - K_active\R_active;
            case 'staggering'%变时间步只修改了staggering
                ElementType = 'P12D';
                Ru = Assemble_half_cohesive(Sol_u,Sol_d,dev_ep_LastLoadStep, IValpha_LastLoadStep, Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                    constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                Rd = Assemble_half_cohesive(Sol_u, Sol_d,dev_ep_LastLoadStep, IValpha_LastLoadStep, Coord, IEN, LM_u, LM_d, elementType, ...   %Rd is the global residual vector for d
                    constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'d');
                Ru_active = Ru(active_indices_u);
                Rd_active = Rd(active_indices_d);
                f_0_t = zeros(length(Ru_active),1);
                conver_vector(iter_no) = norm(Ru_active);
                fprintf('L2-norm of Ru = %g, L2-norm of Rd = %g\n',norm(Ru_active),norm(Rd_active));
                if norm(Ru_active) <= eigenstrain_max * radius * R_tol * sqrt(length(Ru_active)) && norm(Rd_active) <= R_tol * sqrt(length(Rd_active))
                    converged = true;
                    [~,~,dev_ep_new,ep_new_vec,IValpha,crack_RR,plstrain_new_tensor,stress_vector,coord_stressvector] =  Assemble_half_cohesive(Sol_u,Sol_d,...
                        dev_ep_LastLoadStep, IValpha_LastLoadStep,...
                        Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                        constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');%只有收敛了ep等变量才会输出

                    break;
                end
                
                Nt_converged = false;
                for Nt_iter_no1 = 1:max_iter1%u的迭代
                 
                    if Nt_iter_no1 == 1 %删掉
                       [~,Ku,dev_ep_new,ep_new_vec,IValpha,crack_RR,plstrain_new_tensor,stress_vector,coord_stressvector] =  Assemble_half_cohesive(Sol_u,Sol_d, ...
                        dev_ep_LastLoadStep, IValpha_LastLoadStep,...
                        Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                        constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                        Ku_active = Ku(active_indices_u, active_indices_u);
                        Sol_u(active_indices_u) = Sol_u(active_indices_u) - Ku_active\Ru_active;
                        sum_crack_RR = sum(crack_RR);
                        Ru = Assemble_half_cohesive(Sol_u,Sol_d, ...
                            dev_ep_LastLoadStep, IValpha_LastLoadStep,  ...
                            Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                            constitutive, body_force, traction, load_pre,step_no, nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                        Ru_active = Ru(active_indices_u);
                        Ruactive_vector(Nt_iter_no1) = norm(Ru_active);
                        fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
                        if norm(Ru_active) <= eigenstrain_max * radius * R_tol * sqrt(length(Ru_active))
                            Nt_converged = true; break;
                        end
                    else
                        u_c = eigenstrain_max * radius;
                        t_c = deltau;
                        if fixed ==1
                            [~,Ku,dev_ep_new,ep_new_vec,IValpha,crack_RR,plstrain_new_tensor,stress_vector,coord_stressvector] =  Assemble_half_cohesive(Sol_u,Sol_d, ...
                            dev_ep_LastLoadStep, IValpha_LastLoadStep,...
                            Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                            constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                            Ku_active = Ku(active_indices_u, active_indices_u);
                            Sol_u(active_indices_u) = Sol_u(active_indices_u) - Ku_active\Ru_active;
%                             sum_crack_RR = sum(crack_RR);
                            Ru = Assemble_half_cohesive(Sol_u,Sol_d, ...
                                dev_ep_LastLoadStep, IValpha_LastLoadStep,  ...
                                Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                                constitutive, body_force, traction, load_pre,step_no, nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                            Ru_active = Ru(active_indices_u);
                            Ruactive_vector(Nt_iter_no1) = norm(Ru_active);
                            fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
                            if norm(Ru_active) <= eigenstrain_max * radius * R_tol * sqrt(length(Ru_active))
                                delta_u = Sol_u(active_indices_u) - Sol_u_LastLoadStep(active_indices_u);
%                                 delta_R = sum_crack_RR - sum_crack_RR_old;
                                delta_t = load_pre - load_old;
                                if (norm(delta_u)^2/(length(Sol_u(active_indices_u))*u_c^2) + delta_t^2/t_c^2 - 0.1) <= 0                                
                                    Nt_converged = true; break;
                                else
                                    fixed = 0;
                                end
                            end
                        end
                        if fixed == 0
                            delta_u = Sol_u(active_indices_u) - Sol_u_LastLoadStep(active_indices_u);
%                             delta_R = sum_crack_RR - sum_crack_RR_old;
                            delta_t = load_pre - load_old;
                            [~,Ku,dev_ep_new,ep_new_vec,IValpha,crack_RR,plstrain_new_tensor,stress_vector,coord_stressvector] =  Assemble_half_cohesive(Sol_u,Sol_d, ...
                            dev_ep_LastLoadStep, IValpha_LastLoadStep,...
                            Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                            constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                            Ku_active = Ku(active_indices_u, active_indices_u);
                            Ru_active_arclength = [Ru_active;norm(delta_u)^2/(length(Sol_u(active_indices_u))*u_c^2) + delta_t^2/t_c^2 - 0.1];
                            Ku_active_arclength = [[Ku_active,f_0_t];2*delta_u'/(length(Sol_u(active_indices_u))*u_c^2),2*delta_t/t_c^2];
                            delta_arclength = -Ku_active_arclength\Ru_active_arclength;
                            Sol_u(active_indices_u) = Sol_u(active_indices_u) + delta_arclength(1:length(Sol_u(active_indices_u)));%u
                            load_pre = load_pre + delta_arclength(length(Sol_u(active_indices_u))+1);%t
%                             sum_crack_RR = sum(crack_RR);%RR
                            Ru = Assemble_half_cohesive(Sol_u,Sol_d, ...
                                dev_ep_LastLoadStep, IValpha_LastLoadStep,  ...
                                Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                                constitutive, body_force, traction, load_pre,step_no, nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                            Ru_active = Ru(active_indices_u);
                            Ruactive_vector(Nt_iter_no1) = norm(Ru_active);
                            fprintf('L2-norm of Ru = %g\n',norm(Ru_active));
                            if norm(Ru_active) <= eigenstrain_max * radius * R_tol * sqrt(length(Ru_active))
                                Nt_converged = true; break;
                            end
                        end
                    end
                end
                fprintf('load_present is %f \n',load_pre);
                if ~Nt_converged
                    fprintf('There are too many Ru Newton iterations!');
                    converged = false; break;
                end
                
                
                Rd = Assemble_half_cohesive(Sol_u, Sol_d,...
                    dev_ep_LastLoadStep, IValpha_LastLoadStep, ...
                    Coord, IEN, LM_u, LM_d, elementType, ...
                    constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'d');
                Rd_active = Rd(active_indices_d);
                fprintf('L2-norm of Rd = %g\n',norm(Rd_active));
                if norm(Ru_active) <= eigenstrain_max * radius * R_tol * sqrt(length(Ru_active)) && norm(Rd_active) <= R_tol * sqrt(length(Rd_active))
                    converged = true;
                    [~,~,dev_ep_new,ep_new_vec,IValpha,crack_RR,plstrain_new_tensor,stress_vector,coord_stressvector] =  Assemble_half_cohesive(Sol_u,Sol_d, ...
                        dev_ep_LastLoadStep, IValpha_LastLoadStep,...
                        Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                        constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                    break;
                end
                
                Nt_converged = false;
                for Nt_iter_no2 = 1:max_iter2%d的迭代
                    [~,Kd] = Assemble_half_cohesive(Sol_u, Sol_d, dev_ep_LastLoadStep, IValpha_LastLoadStep, Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                        constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'d');
                    Kd_active = Kd(active_indices_d, active_indices_d);
                    Sol_d(active_indices_d) = Sol_d(active_indices_d) - Kd_active\Rd_active;
                    Rd = Assemble_half_cohesive(Sol_u, Sol_d, dev_ep_LastLoadStep, IValpha_LastLoadStep, Coord, IEN, LM_u, LM_d, elementType, ...   %Rd is the global residual vector for d
                        constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'d');
                    Rd_active = Rd(active_indices_d);
                    %
                    Ru = Assemble_half_cohesive(Sol_u,Sol_d, ...
                        dev_ep_LastLoadStep, IValpha_LastLoadStep,  ...
                        Coord, IEN, LM_u, LM_d, elementType, ...   %Ru is the global residual vector for u
                        constitutive, body_force, traction, load_pre,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,'u');
                    Ru_active = Ru(active_indices_u);
                    Ruactive_vector(Nt_iter_no1) = norm(Ru_active);
                    fprintf('L2-norm of Rd = %g, L2-norm of Ru = %g\n',norm(Rd_active),norm(Ru_active));
                    if norm(Rd_active) <= R_tol * sqrt(length(Rd_active))
                        Nt_converged = true; break;
                    end
                end
                if ~Nt_converged
                    fprintf('There are too many Rd Newton iterations!\n');
                    converged = false; break;
                end
        end
    end

    if (converged)
        converge_laststep = true;
        %update the history variables
        Psi_plus_old = Psi_plus_rec;
        Psi_plus_old_d = Psi_plus_rec_d;
        alpha_bar_old = alpha_bar_new;
        pl_energy_all_old = pl_energy_all;
        dev_ep_LastLoadStep = dev_ep_new;
        ep_vector = ep_new_vec;
        IValpha_LastLoadStep = IValpha;
        %         plStrainvector_LastLoadStep =  plStrainvector;
        plstrain_old_tensor = plstrain_new_tensor;
        Sol_u_LastLoadStep = Sol_u;
        Sol_d_LastLoadStep = Sol_d;
        stressvec_old = stress_vector;
        coord_stressvector_old = coord_stressvector;
        disp('Converged.and now output the results');
        if mod(load_pre,period) == 0
            H_gauss(:,fix(load_pre/period)+1) = H_gauss_record;
            alpha_bar_gauss(:,fix(load_pre/period)+1) = alpha_bar_gauss_record;
        end
        %% if converged, output the results
%         load test.mat
        F = Ru(BCIndices_u);
        sum(abs(F));
        if strcmp(iterationType,'staggering')
            Sol(u_indices) = Sol_u;
            Sol(d_indices) = Sol_d;
        end
        
        switch problemType
            case 'PhaseField'
                Answer=reshape(Sol, 3, []);
                Displace=Answer(1:2,:);
                D=Answer(3,:);
            case 'Poisson'
                D=Sol';
        end
        
        
        
        if (abs(Sol_d(526)-1) < 1e-1 && r1 == 0) || (Sol_d(526)-1 > 1e-5 && r1 == 0)
            Loadstep_recorder(1,:) = [step_no fix(load_pre/period)+1];
            r1 = r1+1;
        end


        top=false(2*nNodes,1);
        count = 0;
        for iNode = 1:nNodes
            if Coord(2,iNode) == 0.0000005 
                top(2*iNode)=true;%选取y
                count = count + 1;
            end
        end
        sumRu(step_no,:) = [sum(Ru(top)) load_pre];%displacement control


        % [count_Ru,Ru_sum,Ru_average,Coord_vector,RuCoordR_vector] = RuOfEachNode(Coord,Ru,nNodes);
        % fprintf('count_Ru = %g\n',count_Ru);
        % fprintf('Ru_sum = %g\n',Ru_sum);
        % fprintf('Ru_average = %g\n',Ru_average);
        % fprintf('load_present is %f \n',load_pre);


        
        fprintf('sumRu = %g\n',sum(Ru(top)));
        %输出R(u)的加和
        fprintf('sum of all R(u) = %g\n',sum(Ru));
        %输出R(u)正值的加和
        sum_plus = 0;
        for k = 1 : size(Ru,1)
            if Ru(k) > 0
                sum_plus = sum_plus + Ru(k);
            end
        end
        fprintf('sum_plus = %g\n',sum_plus);
        
      
        
        X = [Coord(1,IEN(1,:)); Coord(1,IEN(2,:)); Coord(1,IEN(3,:))];
        Y = [Coord(2,IEN(1,:)); Coord(2,IEN(2,:)); Coord(2,IEN(3,:))];
        D_elem = [D(IEN(1,:)); D(IEN(2,:)); D(IEN(3,:))];
      

        period_time=2;
        cycle=fix(load_pre/period_time)+1;
       

        ep_elem = [ep_new_vec(1,:); ep_new_vec(2,:); ep_new_vec(3,:)];    %%plStrainvector in yy direction

        
        
        c_length(step_no) = sum(crack_RR);
        if mod(load_pre,period) == 0
            c_length_recorder(fix(load_pre/period)+1,:) = [fix(load_pre/period)+1 sum(crack_RR) sum(sum(Psi_plus_old)) sum(sum(Psi_plus_old_d)) sum(sum(driving_force))];
        end
        

        fprintf('c_length(%g) = %g under t= %g\n',step_no,c_length(step_no),load_pre);


  


        sum_alpha(step_no) = sum(sum(alpha_bar_new));
        fprintf('sum_alpha(%g) = %g under t= %g\n',step_no,sum_alpha(step_no),load_pre);
        norm_alpha(step_no) = norm(alpha_bar_new);
        fprintf('norm_alpha(%g) = %g under t= %g\n',step_no,norm_alpha(step_no),load_pre);
        sum_psi_plus_rec(step_no) = sum(sum(Psi_plus_rec));
        fprintf('sum_psi_plus_rec(%g) = %g under t= %g\n',step_no,sum_psi_plus_rec(step_no),load_pre);

        loadpre_vector(step_no,1) = load_pre;
        %alpha_matrix(:,:,step_no) = alpha_bar_new;
        %d_matrix(:,:,step_no) = D_elem;

        figure(4)
        patch(X,Y,D_elem);
        axis equal;
        axis off;
        colorbar;
        caxis([0 1]);
        title(sprintf('The d under loading step %d ', step_no));


        figure(9)
        norm_alpha_vector(step_no,2) = load_pre;
        norm_alpha_vector(step_no,1) = norm_alpha(step_no);
        scatter(norm_alpha_vector(:,2),norm_alpha_vector(:,1),'filled');
        title(sprintf('The norm of alpha under t= %.6f ', load_pre));
        xlabel('t') ;
        ylabel('norm of alpha') ;
        grid on;



        psi_plus_vector(step_no,1) = sum(sum(Psi_plus_old));
        psi_plus_vector(step_no,2) = load_pre;

        
%         subplot(2,2,3);
        alpha_bar_vector(step_no,1) =  sum(sum(alpha_bar_old));
        alpha_bar_vector(step_no,2) = load_pre;
%         scatter(alpha_bar_vector(:,2),alpha_bar_vector(:,1));
%         title(sprintf('The alpha_bar_sum under loading step %d ', step_no));
%         xlabel('defomation') ;
%         ylabel('alpha_bar_sum') ;
        



        pic_num=pic_num+1;
        

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
        if load_pre<5
        filename = 'test.mat';
        save(filename);
        end
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
end
% save ('Ru_average.mat', 'Ru_record');




