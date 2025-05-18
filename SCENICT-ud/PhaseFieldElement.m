function [PI_ud, crack_r,r_e, k_e,u_DoFs,d_DoFs,Stress_quads,coord_Stress_quads] = PhaseFieldElement(localCoord, localSol, LOCALSol,   ...
    ElementType, constitutive, body_force, traction, current_load_step,step_no, ielem,local_gc)
%	Calculates the element residual vector and the element tangent stiffness matrix
%   r_e = int (sigma(Sol) * epsilon(Na)) - int (Body Force * Na)
%       - int (traction * Na) on Neumann boundary
%   k_e = int (epsilon(Na) * elasticity_matrix(Sol) * epsilon(Na))
% Inputs:
%   localCoord: nDim * nNodesElement
%   LOCALSol: (ux,uy,d)
%       The order of localSol is always (u_x, u_y, (u_z), d) for node 1, then
%       (u_x, u_y, (u_z), d) for node 2, etc.
%   ElementType: 'P12D', etc.
% Output:
%   r_e: (nDoF*nNodesElement) * 1
%   k_e: (nDoF*nNodesElement) * (nDoF*nNodesElement)
% Intermediates:
%   StrainVector, StressVector: 1*1 for 1D, 3*1 for 2D, 6*1d for 3D
%       2D indices: 3=xy
%       3D indices: 4=yz, 5=xz, 6=xy
%   BodyForces: nDim*1
%   detJ: 1 * nQuad
%   Na, dNa_dx, dNa_dy, dNa_dz: nNodesElement * nQuad
%   B: {1,3,6} * (nNodesElement * nDim)
if (nargin < 15)
    switch ElementType
        case 'P12D'
            nQuad = 3; nQuadBdry = 2;
        case 'Q12D'
            nQuad = 4; nQuadBdry = 2;
    end
end
[xi, w] = GetQuadratureRule(ElementType, nQuad);
nDim = size(localCoord, 1);%localCoord2*3
nDoF = nDim + 1;
nNodesElement = size(localCoord, 2);%3
global Psi_plus_rec Psi_plus_old 
% global alpha_bar_old alpha_bar_new
% global stressvector_gauss_y stressvector_gauss_x stressvector_gauss_xy Hydrostatic_stress %strainvector_theta_gauss %µ¯ÐÔÓ¦±ä
% global R_gauss x_gauss y_gauss strainvector_gauss_x strainvector_gauss_y 
%  global lambda_max
% % global coord_stressvector_guass_R coord_stressvector_guass_theta coord_stressvector_guass_psi
% global eigen_strain_matrix
% global strain_x strain_y strain_xy  strain_z
% % global f_alpha_bar
% % global pl_energy_all
% global H_gauss_record global_lambda_max
% global temperature

%u&d for LOCALSol, LOCALSol is used for assemble, K is just for u&d
d_DoFs = nDoF * (1:nNodesElement)';     %3  6   9
u_DoFs = true(length(LOCALSol), 1);
u_DoFs(d_DoFs) = false;     %110 110 110
T_DoFs = (nDoF+1) * (1:nNodesElement)';
% T_DoFs = d_DoFs+1;
% %used for sol, sol includes ux,uy,d,T
% local_d_indices = nDim+1 : nDoF : (nDoF * nNodesElement);   %3  7   11
% local_T_indices = nDoF : nDoF : (nDoF * nNodesElement);     %4  8   12
% local_u_indices = setdiff(1:(nDoF * nNodesElement), local_d_indices);
% local_u_indices = setdiff(local_u_indices, local_T_indices);    %12 56 910
% uy_on_gauss = zeros(nQuad,1);
% ep_new_local_vec = zeros(nQuad,1);
 Stress_quads = zeros(3 *  nQuad,1);
%Stress_quads = zeros(4 *  nQuad,1);%sys12*1
coord_Stress_quads = zeros(3 *  nQuad,1);
% EleStressvector = zeros(nNodesElement*nDoF,1);
% EleStrainvector = zeros(nNodesElement*nDoF,1);

switch nDim
    case 1
        [detJ, Na, dNa_dx] = QuadShape(ElementType, localCoord, xi);
    case 2
        [detJ, Na, dNa_dx, dNa_dy] = QuadShape(ElementType, localCoord, xi);
    case 3
        [detJ, Na, dNa_dx, dNa_dy, dNa_dz] = QuadShape(ElementType, localCoord, xi);
end

r_e_u = zeros(nDim*nNodesElement, 1);
r_e_d = zeros(nNodesElement, 1);
crack_length = zeros(nNodesElement, 1);
crack_length_sum = zeros(nNodesElement, 1);
k_e_uu = zeros(length(r_e_u));
k_e_dd = zeros(length(r_e_d));
k_e_ud = zeros(length(r_e_u), length(r_e_d));
zero_filler = zeros(1, length(dNa_dx));
PI_ud = 0;
for iQuad = 1:nQuad
    x_quad = localCoord(1,:) * Na(:,iQuad); %x coordinate of gauss point
%     integral_factor = 2 * pi * r_quad;
    y_quad = localCoord(2,:) * Na(:,iQuad); %y coordinate of gauss point
%     R_quad = sqrt(x_quad^2+y_quad^2);  
    T = localSol(T_DoFs)' * Na(:,iQuad);    %Temperature of gauss point

%     coord_psi=asin(x_quad/R_quad);
%     coord_theta=0;
    
%     coord_psi=atan(abs(r_quad)/z_quad);
%     coord_psi_matrix(iQuad,ielem)=coord_psi*180/pi;
    switch nDim
        case 1
            B = dNa_dx(:,iQuad);
            Bd = B;
        case 2
            B = [reshape([dNa_dx(:,iQuad)'; zero_filler], 1, [])
                reshape([zero_filler; dNa_dy(:,iQuad)'], 1, [])
                reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'], 1, [])];
            Bd = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'];  %2*3
            
%             B_sys = [reshape([dNa_dx(:,iQuad)'; zero_filler], 1, [])
%                 reshape([Na(:,iQuad)'./x_quad; zero_filler], 1, [])
%                 reshape([zero_filler; dNa_dy(:,iQuad)'], 1, [])
%                 reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'], 1, [])];
%             Bd_sys = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'];
            
        case 3
            B = [reshape([dNa_dx(:,iQuad)'; zero_filler; zero_filler], 1, [])
                reshape([zero_filler; dNa_dy(:,iQuad)'; zero_filler], 1, [])
                reshape([zero_filler; zero_filler; dNa_dz(:,iQuad)'], 1, [])
                reshape([zero_filler; dNa_dz(:,iQuad)'; dNa_dy(:,iQuad)'], 1, [])
                reshape([dNa_dz(:,iQuad)'; zero_filler; dNa_dx(:,iQuad)'], 1, [])
                reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'; zero_filler], 1, [])];
            Bd = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'; dNa_dz(:,iQuad)'];
    end
    %epsilony_AT2 = sqrt(gc / (3 * ell * E));
    %alphat = 1/2 * E * epsilony_AT2 * epsilony_AT2;%N/mm^2
%     alphat =40;
%     ka = 1/3;
    
    %     StrainVector = B * localSol(u_DoFs);
    %     StrainVector = B_sys * localSol(u_DoFs);
    [ei_strain] = eigen_strain(T);
%     eigen_strain_matrix(iQuad,ielem)= ei_strain;
%     if ei_strain~=0
%         i=1;
%     end
    StrainVector = B * LOCALSol(u_DoFs) - ei_strain*[1;1;0];
%eigen_strain(R_quad,current_load_step)
    %elastic strain      total strain  3*6 6*1                eigen strain
    %StrainVector_t = -ei_strain_t;
    %fprintf('StrainVector(1)= %g\n',StrainVector(1));
    strain = StrainVector + ei_strain*[1;1;0];  %total strain
%     strain_x(iQuad,ielem)=strain(1);
%     strain_y(iQuad,ielem)=strain(2);
%     strain_xy(iQuad,ielem)=strain(3);

%     strain_rz(iQuad,ielem)=strain(4);
%     localSol_cur = LOCALSol(u_DoFs);
%     quad_gc = local_gc *  Na(:,iQuad);
%     ep_old = ep_old_vector(iQuad);

    
    d = LOCALSol(d_DoFs)' * Na(:,iQuad);
    grad_d = Bd * LOCALSol(d_DoFs); %2*1
    %modelB
%     [StressVector, D, Psi_plus, sigma_plus, ~, gc, ell, ~, local_dev_ep_new((iQuad-1) * 6 + 1 : iQuad * 6), localIValpha_new((iQuad-1) * 6 + 1 : iQuad * 6),~,ep_new,plstrain_new_vector((iQuad-1) * 6 + 1 : iQuad * 6),pl_strain_energy] =...
%         constitutive(StrainVector,d,local_dev_ep_LastLoadStep((iQuad-1) * 6 + 1 : iQuad * 6), localIValpha_LastLoad((iQuad-1) * 6 + 1 : iQuad * 6),quad_gc,ep_old,plstrain_old_vector((iQuad-1) * 6 + 1 : iQuad * 6));
       %----modelC----%
[StressVector, D, Psi, Psi_plus, sigma_plus, gc,ell] =...
        constitutive(StrainVector,T, d, y_quad ,'PlaneStress');
% local_dev_ep_new = local_dev_ep_LastLoadStep;
% localIValpha_new = localIValpha_LastLoad;
% plstrain_new_vector = plstrain_old_vector;
           %----modelC----%

gd = (1-d)^2 + 1e-10;
% ep_new = 0;
%     ep_new_local_vec(iQuad) = ep_new;
    Stress_quads(iQuad*3-2:iQuad*3) = StressVector;

    StressMatrix = [StressVector(1),StressVector(3)/2,0;
        StressVector(3)/2, StressVector(2),0;
        0 ,0 , 0]; 
    lambda = eig(StressMatrix);
%     global_lambda_max(iQuad,ielem) = max(max(lambda));
%     if step_no == 1
%         if d<exp(-2)
%             lambda_max(iQuad,ielem) = max(max(lambda));
%         else 
%             lambda_max(iQuad,ielem) = 0;
%         end
%     else 
%         lambda_max(iQuad,ielem) = max(max(lambda));
%     end
%     
%     stressvector_gauss_x(iQuad,ielem) = StressVector(1);
%     stressvector_gauss_y(iQuad,ielem) = StressVector(2);
%     stressvector_gauss_xy(iQuad,ielem) = StressVector(3);
%     Hydrostatic_stress(iQuad,ielem) = 1/3*(StressVector(1)+StressVector(2)+StressVector(3));
% 
% %     strainvector_theta_gauss(iQuad,ielem) = StrainVector(2);
% %     R_gauss(iQuad,ielem) = R_quad;
%     x_gauss(iQuad,ielem) = x_quad;
%     y_gauss(iQuad,ielem) = y_quad;
%     strainvector_gauss_x(iQuad,ielem) = StrainVector(1);
%     strainvector_gauss_y(iQuad,ielem) = StrainVector(2);
%     beta_1 = [sin(coord_psi)*cos(coord_theta),sin(coord_psi)*sin(coord_theta),cos(coord_psi)
%         -sin(coord_theta),cos(coord_theta),0
%         cos(coord_psi)*cos(coord_theta),sin(coord_psi)*sin(coord_theta),-cos(coord_psi)];  
%     beta_2 = [cos(coord_theta),sin(coord_theta),0
%         -sin(coord_theta),cos(coord_theta),0
%         0,0,1];
%     beta_2_inv=inv(beta_2);
%     beta_2_t_inv=inv(beta_2.');
%     StressMatrix = [StressVector(1),StressVector(3)/2,0;
%         StressVector(3)/2, StressVector(2),0;
%         0 ,0 , 0]; 
%     sigma_d = beta_2_inv*StressMatrix*beta_2_t_inv;     %Radial stress
%     sigma_s = beta_1*sigma_d*beta_1.';  %Circumferential stress
%     coord_stressvector = sigma_s;%3*3
%     lambda = eig(StressMatrix);
%     lambda_1(iQuad,ielem) = lambda(1);
%     lambda_2(iQuad,ielem) = lambda(2);
%     lambda_3(iQuad,ielem) = lambda(3);
%     lambda_max(iQuad,ielem) = max(lambda(1),max(lambda(2),lambda(3)));

    %Spherical coordinates
%     coord_stressvector=[sin(coord_psi),0,cos(coord_psi);0,1,0;cos(coord_psi),0,-sin(coord_psi)]*[StressVector(1),StressVector(2),StressVector(3)].';
%     coord_Stress_quads(iQuad*3-2:iQuad*3)= [coord_stressvector(1);coord_stressvector(5);coord_stressvector(9)];
%     coord_stressvector_guass_R(iQuad,ielem)=coord_stressvector(1);
%     coord_stressvector_guass_theta(iQuad,ielem)=coord_stressvector(5);
%     coord_stressvector_guass_psi(iQuad,ielem)=coord_stressvector(9);

    %To prevent crack healing, psi is the energy driving crack growth, so this value can only increase but not decrease
%     Psi_th = 180e6^2/(2*370e9); %material strength^2/2E
    Psi_plus_rec(iQuad,ielem) = max(Psi_plus_old(iQuad,ielem),Psi_plus); %record for history
%     Pure_Psi_plus(iQuad,ielem)= Psi_plus;
    Psi_plus_new =  Psi_plus_rec(iQuad,ielem);
%     Psi_plus_rec_d(iQuad,ielem) =  gd * Psi_plus;
%     H_gauss_record(iQuad + (ielem - 1) * nQuad) = gd * Psi_plus;
%     driving_force(iQuad,ielem) = 2 * (1-d)* Psi_plus_new;
    
%     f_alpha_bar(iQuad,ielem) = 1;
%     %history variable£ºFatigue degradation function
%     %alpha_bar_new(iQuad,ielem) = alpha_bar_old(iQuad,ielem) + (Psi_plus_new - Psi_plus_old(iQuad,ielem) + abs(etat * ep_new))*heaviside(Psi_plus_new - Psi_plus_old(iQuad,ielem)+abs(etat * ep_new));
%     alpha_bar_new(iQuad,ielem) = alpha_bar_old(iQuad,ielem) + (Psi_plus_rec_d(iQuad,ielem) - Psi_plus_old_d(iQuad,ielem))*heaviside(Psi_plus_rec_d(iQuad,ielem) - Psi_plus_old_d(iQuad,ielem));
%     %         alpha_bar_new(iQuad,ielem) = alpha_bar_old(iQuad,ielem) + (Psi_plus_rec_d(iQuad,ielem) - Psi_plus_old_d(iQuad,ielem))*heaviside(Psi_plus_rec_d(iQuad,ielem) - Psi_plus_old_d(iQuad,ielem));
%     alpha_bar_gauss_record(iQuad + (ielem - 1) * nQuad) = alpha_bar_new(iQuad,ielem);
% 
%     if  alpha_bar_new(iQuad,ielem) < alphat 
%         f_alpha_bar(iQuad,ielem) = 1;
%     end
%     if alpha_bar_new(iQuad,ielem) < alphat * 10^(1/ka) && alpha_bar_new(iQuad,ielem) > alphat
%         f_alpha_bar(iQuad,ielem) = (1 - ka * log10(alpha_bar_new(iQuad,ielem)/alphat))^2;
%     end
%     if  alpha_bar_new(iQuad,ielem) > alphat * 10^(1/ka)
%         f_alpha_bar(iQuad,ielem) = 0;
%     end
    %     if quad_gc == 1000
    %         f_alpha_bar(iQuad,ielem) = 1;
    %     end
%     f_alpha_bar(iQuad,ielem) = 1;
      BodyForces = body_force(localCoord * Na(:,iQuad));  %bodyforce of gauss point
%     BodyForces = current_load_step * body_force(localCoord * Na(:,iQuad));  %bodyforce of gauss point,can grow over time
    r_e_u = r_e_u + (B' * StressVector - reshape(BodyForces * Na(:,iQuad)', [], 1)) ...
        * w(iQuad) * detJ(iQuad);
%         * integral_factor;%B_sys'6*4 
    %     r_e_d = r_e_d + (- 2 * (1-d) * Na(:,iQuad) * (Psi_plus_new + 0.1 * pl_energy_all(iQuad,ielem)) + f_alpha_bar(iQuad,ielem) * gc * (d * Na(:,iQuad) / ell ...
    %         + ell * Bd' * grad_d)) ...
    %         * w(iQuad) * detJ(iQuad);
    r_e_d = r_e_d + (- 2 * (1-d) * Na(:,iQuad) * (Psi_plus_new) + gc * (d * Na(:,iQuad) / ell ...
        + ell * Bd' * grad_d)) ...
        * w(iQuad) * detJ(iQuad);
%fatigue
%     r_e_d = r_e_d + (- 2 * (1-d) * Na(:,iQuad) * (Psi_plus_new + 0 * pl_energy_all(iQuad,ielem)) + (f_alpha_bar(iQuad,ielem)) * gc * (d * Na(:,iQuad) / ell ...
%         + ell * Bd' * grad_d)) ...
%         * w(iQuad) * detJ(iQuad);
%         * integral_factor;
    crack_length(iQuad) = ((d * d / (2 * ell) ...
        + ell / 2 * (grad_d' * grad_d))) ...
        * w(iQuad) * detJ(iQuad);
%         * integral_factor;


    k_e_uu = k_e_uu + B' * D * B* w(iQuad) * detJ(iQuad);
%     k_e_dd = k_e_dd + (Na(:,iQuad) * (2 * (Psi_plus_new +   0 * pl_energy_all(iQuad,ielem)) + (f_alpha_bar(iQuad,ielem)) * gc / ell) * Na(:,iQuad)' ...
%         + Bd' * (f_alpha_bar(iQuad,ielem)) * gc * ell * Bd) ...
%         * w(iQuad) * detJ(iQuad);
%         * integral_factor;
            k_e_dd = k_e_dd + (Na(:,iQuad) * (2 * (Psi_plus_new) + gc / ell) * Na(:,iQuad)' ...
                + Bd' * gc * ell * Bd) ...
                * w(iQuad) * detJ(iQuad);
    k_e_ud = k_e_ud - 2 * B' * (1-d) * sigma_plus * Na(:,iQuad)' ...
        * w(iQuad) * detJ(iQuad);
    PI_ud = PI_ud + (Psi + gc * (d^2+ell^2*abs(grad_d'*grad_d))/(2*ell))* w(iQuad) * detJ(iQuad);

%         * integral_factor;
end

switch ElementType
    case 'P12D'
        nFaces = 3;
        FaceNodes = [1, 2, 3
            2, 3, 1];
        FaceElementType = 'P11D';
    case 'Q12D'
        nFaces = 4;
        FaceNodes = [1, 2, 3, 4
            2, 3, 4, 1];
        FaceElementType = 'P11D';
    case 'P13D'
        nFaces = 4;
        FaceNodes = [1, 1, 1, 2
            3, 2, 4, 3
            2, 4, 3, 4];
end

% Compute the contribution of the traction to the residual.
[xi, w] = GetQuadratureRule(FaceElementType, nQuadBdry);
for iFace = 1:nFaces%3
    ThisFaceNodes = FaceNodes(:, iFace);
    FaceCoord = localCoord(:, ThisFaceNodes);
    element_uDoFs = reshape(repmat(ThisFaceNodes'-1, nDim, 1) * nDim + repmat((1:nDim)', 1, length(ThisFaceNodes)), [], 1);
    [detJ, Na] = FaceQuadShape(FaceElementType, FaceCoord, xi);
    for iQuadBdry = 1:nQuadBdry%2
        TractionThisPoint = traction(FaceCoord * Na(:, iQuadBdry),step_no) ;%* integral_factor;
        r_e_u(element_uDoFs) = r_e_u(element_uDoFs) - reshape(TractionThisPoint * Na(:, iQuadBdry)', [], 1) ...
            * w(iQuadBdry) * detJ(iQuadBdry);
    end
end

r_e = zeros(length(LOCALSol), 1);
crack_r = zeros(length(LOCALSol), 1);
r_e(u_DoFs) = r_e_u;
r_e(d_DoFs) = r_e_d;
crack_r(d_DoFs) =  crack_length;

k_e = zeros(length(LOCALSol), length(LOCALSol));
k_e(u_DoFs, u_DoFs) = k_e_uu;
k_e(d_DoFs, d_DoFs) = k_e_dd;
k_e(u_DoFs, d_DoFs) = k_e_ud;
k_e(d_DoFs, u_DoFs) = k_e_ud';