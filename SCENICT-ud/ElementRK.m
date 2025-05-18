function [r_e, k_e] = ElementRK(localCoord, localSol,   ...
    ElementType, constitutive, body_force, traction, localSol_T_record ,step_no, ielem)
switch ElementType
    case 'P12D'
        nQuad = 3; nQuadBdry = 2;
    case 'Q12D'
        nQuad = 4; nQuadBdry = 2;
end
[xi, w] = GetQuadratureRule(ElementType, nQuad);
nDim = size(localCoord, 1);%localCoord2*3
nDoF = nDim + 2;%位移场，相场，温度场
nNodesElement = size(localCoord, 2);%3
global Psi_plus_rec Psi_plus_old 
global deltau T0
T_DoFs = nDoF * (1:nNodesElement)';     %4   8   12
d_DoFs = nDoF * (1:nNodesElement)' - 1; %3   7   11
u_DoFs = true(length(localSol),1);
u_DoFs(d_DoFs) = false;
u_DoFs(T_DoFs) = false;
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
r_e_T = zeros(nNodesElement, 1);
k_e_uu = zeros(length(r_e_u));
k_e_dd = zeros(length(r_e_d));
k_e_TT = zeros(length(r_e_T));
k_e_ud = zeros(length(r_e_u), length(r_e_d));
k_e_du = zeros(length(r_e_d), length(r_e_u));
k_e_uT = zeros(length(r_e_u), length(r_e_T));
k_e_dT = zeros(length(r_e_d), length(r_e_T));
k_e_Td = zeros(length(r_e_T), length(r_e_d));
zero_filler = zeros(1, length(dNa_dx));
for iQuad = 1:nQuad
    x_quad = localCoord(1,:) * Na(:,iQuad); %x coordinate of gauss point
    y_quad = localCoord(2,:) * Na(:,iQuad); %y coordinate of gauss point
    T = localSol(T_DoFs)' * Na(:,iQuad);    %Temperature of gauss point

    switch nDim
        case 1
            B = dNa_dx(:,iQuad);
            Bd = B;
        case 2
            B = [reshape([dNa_dx(:,iQuad)'; zero_filler], 1, [])
                reshape([zero_filler; dNa_dy(:,iQuad)'], 1, [])
                reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'], 1, [])];
            Bd = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'];  %2*3
            BT = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'];  %2*3

        case 3
            B = [reshape([dNa_dx(:,iQuad)'; zero_filler; zero_filler], 1, [])
                reshape([zero_filler; dNa_dy(:,iQuad)'; zero_filler], 1, [])
                reshape([zero_filler; zero_filler; dNa_dz(:,iQuad)'], 1, [])
                reshape([zero_filler; dNa_dz(:,iQuad)'; dNa_dy(:,iQuad)'], 1, [])
                reshape([dNa_dz(:,iQuad)'; zero_filler; dNa_dx(:,iQuad)'], 1, [])
                reshape([dNa_dy(:,iQuad)'; dNa_dx(:,iQuad)'; zero_filler], 1, [])];
            Bd = [dNa_dx(:,iQuad)'; dNa_dy(:,iQuad)'; dNa_dz(:,iQuad)'];
    end

    [ei_strain] = eigen_strain(T);
    StrainVector = B * localSol(u_DoFs) - ei_strain*[1;1;0]; %Elastic strain
%     strain = StrainVector + ei_strain*[1;1;0];  %total strain
    d = localSol(d_DoFs)' * Na(:,iQuad);
    grad_d = Bd * localSol(d_DoFs); %2*1
    grad_T = BT * localSol(T_DoFs); %2*1
    [StressVector, D, Psi, Psi_plus, sigma_plus, gc,ell, psi_plus_T] =...
        constitutive(StrainVector,T, d, y_quad ,'PlaneStress');
    gd = (1-d)^2 + 1e-10;
    T_record = localSol_T_record' * Na(:,iQuad);
    rho = 3980;  %Density
    Cp = 880; %Specific heat capacity
    k0 = 31; %inherent thermal conductivity of the material
    Ea = [7.5e-6;7.5e-6;0];%需要和eigen_strain代码里的alpha一样


    Psi_plus_rec(iQuad,ielem) = max(Psi_plus_old(iQuad,ielem),Psi_plus); %record for history



    Psi_plus_new =  Psi_plus_rec(iQuad,ielem);
    BodyForces = body_force(localCoord * Na(:,iQuad));  %bodyforce of gauss point

    r_e_u = r_e_u + (B' * StressVector - reshape(BodyForces * Na(:,iQuad)', [], 1)) ...
        * w(iQuad) * detJ(iQuad);
    r_e_d = r_e_d + (- 2 * (1-d) * Na(:,iQuad) * (Psi_plus_new) + gc * (d * Na(:,iQuad) / ell ...
        + ell * Bd' * grad_d)) ...
        * w(iQuad) * detJ(iQuad);
    r_e_T = r_e_T + (rho * Cp * ((T - T_record)/deltau) *  Na(:,iQuad) + (1 - d)^2 * k0 * BT' * grad_T)...
        * w(iQuad) * detJ(iQuad);

    k_e_uu = k_e_uu + B' * D * B* w(iQuad) * detJ(iQuad);
    k_e_dd = k_e_dd + (Na(:,iQuad) * (2 * (Psi_plus_new) + gc / ell) * Na(:,iQuad)' ...
        + Bd' * gc * ell * Bd) ...
        * w(iQuad) * detJ(iQuad);
    k_e_TT = k_e_TT + (((rho * Cp)/deltau) * Na(:,iQuad) * Na(:,iQuad)' ...
        + gd * k0 * (BT' * BT))* w(iQuad) * detJ(iQuad);
    k_e_ud = k_e_ud - 2 * B' * (1-d) * sigma_plus * Na(:,iQuad)' ...
        * w(iQuad) * detJ(iQuad);
    if Psi_plus_rec(iQuad,ielem)>Psi_plus_old(iQuad,ielem)

        k_e_du = k_e_du - (2 * B' * (1-d) *  sigma_plus * Na(:,iQuad)')' ...
            * w(iQuad) * detJ(iQuad);
%         k_e_ud = k_e_du';
        k_e_dT = k_e_dT - 2 * (1 - d) * Na(:,iQuad) * psi_plus_T * Na(:,iQuad)'...
            * w(iQuad) * detJ(iQuad);

    end
    
    k_e_uT = k_e_uT - B' * D * Ea * Na(:,iQuad)'* w(iQuad) * detJ(iQuad);
    
    k_e_Td = k_e_Td - 2 * (1 - d) * k0 * BT' * grad_T * Na(:,iQuad)'* w(iQuad) * detJ(iQuad);
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

r_e = zeros(length(localSol), 1);

r_e(u_DoFs) = r_e_u;
r_e(d_DoFs) = r_e_d;
r_e(T_DoFs) = r_e_T;

k_e = zeros(length(localSol), length(localSol));
k_e(u_DoFs, u_DoFs) = k_e_uu;
k_e(d_DoFs, d_DoFs) = k_e_dd;
k_e(T_DoFs, T_DoFs) = k_e_TT;
k_e(u_DoFs, d_DoFs) = k_e_ud;
k_e(d_DoFs, u_DoFs) = k_e_du;
k_e(u_DoFs, T_DoFs) = k_e_uT;
k_e(d_DoFs, T_DoFs) = k_e_dT;
k_e(T_DoFs, d_DoFs) = k_e_Td;

end

