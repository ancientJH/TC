function  [R, K] = Assemble_half_cohesive(Sol_u,Sol_d,Sol_T,Coord, IEN, LM_u, LM_d, LM_T, elementType, ...
    constitutive,body_force,traction, current_load_step,step_no,nQuad,gc_vec,WhichHalf)
% function  [R, K, dev_ep_new,ep_new_vec, IValpha_new, crack_RR,stress_vector,coord_stressvector] = Assemble_half_cohesive(Sol_u,Sol_d, dev_ep_LastLoadStep, IValpha_LastLoadStep, Coord, IEN, LM_u, LM_d, elementType, ...
%     constitutive,body_force,traction, current_load_step,step_no,nQuad,gc_vec,ep_vector,plstrain_old_tensor,WhichHalf)
nNodesElement = size(IEN, 1);
nElements = size(IEN, 2);
nDim = size(Coord,1);
nNodes = size(Coord, 2);
nEquations_u = length(Sol_u);
nEquations_d = length(Sol_d);
nEquations_T = length(Sol_T);
nDoF = nDim+2;
% nNode = size(Coord, 2);
% RuCoord = zeros(nDim,nNode);
% ep_new_vec = zeros(size(IEN));
%plstrain_new_tensor = zeros(6 * nQuad,size(IEN,2));
stress_vector = zeros(3*nQuad,size(IEN,2));
coord_stressvector = zeros(3*nQuad,size(IEN,2));
global Sol_T_record
switch WhichHalf
    case 'u'
        R = zeros(size(Sol_u));     %Ru
        if (nargout >= 2)
            K = spalloc(nEquations_u, nEquations_u, 20*nNodes);
        end

    case 'd'
        R = zeros(size(Sol_d));
        if (nargout >= 2)
            K = spalloc(nEquations_d, nEquations_d, 20*nNodes);
        end

    case 'T'
        R = zeros(size(Sol_T));

        K = spalloc(nEquations_T, nEquations_T, 20*nNodes);


    case 'ud'
        R = zeros(size(Sol_u));
        if (nargout >= 2)
            K = spalloc(nEquations_u, nEquations_d, 20*nNodes);
        end

    case 'du'
        R = zeros(size(Sol_d));
        if (nargout >= 2)
            K = spalloc(nEquations_d, nEquations_u, 20*nNodes);
        end

    case 'uT'
        R = zeros(size(Sol_u));
        if (nargout >= 2)
            K = spalloc(nEquations_u, nEquations_T, 60*nNodes);
        end

    case 'dT'
        R = zeros(size(Sol_d));
        if (nargout >= 2)
            K = spalloc(nEquations_d, nEquations_T, 20*nNodes);
        end

    case 'Td'
        R = zeros(size(Sol_T));
        if (nargout >= 2)
            K = spalloc(nEquations_T, nEquations_d, 20*nNodes);
        end
end

local_d_indices = nDim+1 : nDoF : (nDoF * nNodesElement);
local_T_indices = nDoF : nDoF : (nDoF * nNodesElement);
local_u_indices = setdiff(1:(nDoF * nNodesElement), local_d_indices);
local_u_indices = setdiff(local_u_indices, local_T_indices);


for ielem = 1:nElements
    %fprintf('%g th elements\n',ielem);
    local_gc = gc_vec(IEN(:, ielem));   %local！！！！> one element
    localCoord = Coord(:, IEN(:, ielem));   %coordinate of three nodes of one element
    localSol_u = Sol_u(LM_u(:, ielem));
    localSol_d = Sol_d(LM_d(:, ielem));
    localSol_T = Sol_T(LM_T(:, ielem));
    localSol_T_record = Sol_T_record(LM_T(:, ielem));
    localSol = zeros(length(localSol_u) + length(localSol_d) + length(localSol_T));
    localSol(local_u_indices) = localSol_u;
    localSol(local_d_indices) = localSol_d;
    localSol(local_T_indices) = localSol_T;

    if (ielem <= 220000)

        [r_e, k_e] = ElementRK(localCoord, localSol,   ...
            elementType, constitutive, body_force, traction, localSol_T_record ,step_no, ielem);
        r_e_u = r_e(local_u_indices);
        r_e_d = r_e(local_d_indices);
        r_e_T = r_e(local_T_indices);
        k_e_uu = k_e(local_u_indices, local_u_indices);
        k_e_dd = k_e(local_d_indices, local_d_indices);
        k_e_TT = k_e(local_T_indices, local_T_indices);
        k_e_ud = k_e(local_u_indices, local_d_indices);
        k_e_du = k_e(local_d_indices, local_u_indices);
        k_e_uT = k_e(local_u_indices, local_T_indices);
        k_e_dT = k_e(local_d_indices, local_T_indices);
        k_e_Td = k_e(local_T_indices, local_d_indices);
    else
        [r_e, k_e] = PhaseFieldcohesiveElement(localCoord, LOCALSol, localSol_u, LOCALSol_d, elementType, ...
            constitutive, body_force, traction, current_load_step);


    end
    switch WhichHalf
        case 'u'
            R(LM_u(:, ielem)) = R(LM_u(:, ielem)) + r_e_u;
            if (nargout >= 2)
                K(LM_u(:, ielem), LM_u(:, ielem)) = K(LM_u(:, ielem), LM_u(:, ielem)) + k_e_uu;
            end

        case 'd'
            R(LM_d(:, ielem)) = R(LM_d(:, ielem)) + r_e_d;
            if (nargout >= 2)
                K(LM_d(:, ielem), LM_d(:, ielem)) = K(LM_d(:, ielem), LM_d(:, ielem)) + k_e_dd;
            end

        case 'T'
            R(LM_T(:, ielem)) = R(LM_T(:, ielem)) + r_e_T;

            K(LM_T(:, ielem), LM_T(:, ielem)) = K(LM_T(:, ielem), LM_T(:, ielem)) + k_e_TT;


        case 'ud'
            R(LM_u(:, ielem)) = R(LM_u(:, ielem)) + r_e_u;
            K(LM_u(:, ielem), LM_d(:, ielem)) = K(LM_u(:, ielem), LM_d(:, ielem)) + k_e_ud;
        case 'du'
            R(LM_d(:, ielem)) = R(LM_d(:, ielem)) + r_e_d;
            K(LM_d(:, ielem), LM_u(:, ielem)) = K(LM_d(:, ielem), LM_u(:, ielem)) + k_e_du;
        case 'uT'
            R(LM_u(:, ielem)) = R(LM_u(:, ielem)) + r_e_u;
            K(LM_u(:, ielem), LM_T(:, ielem)) = K(LM_u(:, ielem), LM_T(:, ielem)) + k_e_uT;
        case 'dT'
            R(LM_d(:, ielem)) = R(LM_d(:, ielem)) + r_e_d;

            K(LM_d(:, ielem), LM_T(:, ielem)) = K(LM_d(:, ielem), LM_T(:, ielem)) + k_e_dT;
        case 'Td'
            R(LM_T(:, ielem)) = R(LM_d(:, ielem)) + r_e_T;
            if (nargout >= 2)
                K(LM_T(:, ielem), LM_d(:, ielem)) = K(LM_T(:, ielem), LM_d(:, ielem)) + k_e_Td;
            end

    end
end



