function [StressVector, D2d,Psi, Psi_plus, sigma_plus, gc, ell, psi_plus_T] = ModleC(StrainVector,T, d,y_quad,Whichhalf)
E = 370e9;
nu = 0.3;
gc = 42.47;
ell = 7.5e-5;% 采用Wu老师文章，sigmac=180Mpa,所以网格可以采取1.5e-5
alpha = 7.5e-6;
if strcmp(Whichhalf,'PlaneStress')
    E=E*(1+2*nu)/(1+nu)^2;
    nu=nu/(1+nu);
end
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));
gd = (1-d)^2 + 1e-8;
gd_diff = - 2 * (1-d);
angle_bracket_Plus = @(a) (a+abs(a))/2;
angle_bracket_Minus = @(a) (a-abs(a))/2;
switch size(StrainVector, 1)
    case 3
        nDim = 2;
        StrainTensor = [StrainVector(1), StrainVector(3)/2, 0
            StrainVector(3)/2, StrainVector(2), 0
            0, 0, 0];
    case 6
        nDim = 3;
        StrainTensor = [StrainVector(1), StrainVector(6)/2, StrainVector(5)/2
            StrainVector(6)/2, StrainVector(2), StrainVector(4)/2
            StrainVector(5)/2, StrainVector(4)/2, StrainVector(3)];
end

[principal_directions, principal_strains] = eig(StrainTensor);
%注意，principal_directions的第一列才是向量n1，n12代表向量a1的第二个分量
%和组内公众号相互呼应
n11 = principal_directions(1,1);
n12 = principal_directions(2,1);
n13 = principal_directions(3,1);
n21 = principal_directions(1,2);
n22 = principal_directions(2,2);
n23 = principal_directions(3,2);
n31 = principal_directions(1,3);
n32 = principal_directions(2,3);
n33 = principal_directions(3,3);
%A_eps代表应变的转置矩阵，可以参考公众号文章https://mp.weixin.qq.com/s/rORT3N3B4z2PK-4bEQyASQ
%或者组内implementation details for miehe....pdf
A_eps = [n11^2 n12^2 n13^2 n12*n13 n11*n13 n11*n12; 
         n21^2 n22^2 n23^2 n22*n23 n21*n23 n21*n22; 
         n31^2 n32^2 n33^2 n32*n33 n31*n33 n31*n32; 
         2*n21*n31 2*n22*n32 2*n23*n33 n22*n33+n23*n32 n21*n33+n23*n31 n21*n32+n22*n31; 
         2*n11*n31 2*n12*n32 2*n13*n33 n12*n33+n13*n32 n11*n33+n13*n31 n11*n32+n12*n31; 
         2*n11*n21 2*n12*n22 2*n13*n23 n12*n23+n13*n22 n11*n23+n13*n21 n11*n22+n12*n21];
principal_strains = diag(principal_strains);
%4.eig函数求出来的特征方向，特征值均是以矩阵形式输出的？ ANS:yes，所以我们接下来
%在表示特征值的时候，继续用了diag函数把对角线元素取出来。
Psi_plus = lambda/2 * angle_bracket_Plus(trace(StrainTensor))^2 + ...
    mu * sum(angle_bracket_Plus(principal_strains).^2);  %此处利用sum函数求和非常棒
Dn_Plus = [lambda*heaviside(trace(StrainTensor))+2*mu*heaviside(principal_strains(1)) ...
    lambda*heaviside(trace(StrainTensor)) lambda*heaviside(trace(StrainTensor)) 0 0 0; 
    lambda*heaviside(trace(StrainTensor)) lambda*heaviside(trace(StrainTensor))+2*mu*heaviside(principal_strains(2)) ...
    lambda*heaviside(trace(StrainTensor)) 0 0 0; 
    lambda*heaviside(trace(StrainTensor)) lambda*heaviside(trace(StrainTensor)) ...
    lambda*heaviside(trace(StrainTensor))+2*mu*heaviside(principal_strains(3)) 0 0 0 ; 
    0 0 0 mu*Hab_Plus(principal_strains(2),principal_strains(3)) 0 0; 
    0 0 0 0 mu*Hab_Plus(principal_strains(1),principal_strains(3)) 0; 
    0 0 0 0 0 mu*Hab_Plus(principal_strains(1),principal_strains(2))];
Dn_Minus = [lambda*heaviside(-trace(StrainTensor))+2*mu*heaviside(-principal_strains(1)) ...
    lambda*heaviside(-trace(StrainTensor)) lambda*heaviside(-trace(StrainTensor)) 0 0 0; ...
    lambda*heaviside(-trace(StrainTensor)) lambda*heaviside(-trace(StrainTensor))+2*mu*heaviside(-principal_strains(2)) ...
    lambda*heaviside(-trace(StrainTensor)) 0 0 0; ...
    lambda*heaviside(-trace(StrainTensor)) lambda*heaviside(-trace(StrainTensor)) ...
    lambda*heaviside(-trace(StrainTensor))+2*mu*heaviside(-principal_strains(3)) 0 0 0 ; ...
    0 0 0 mu*Hab_Minus(principal_strains(2),principal_strains(3)) 0 0; ...
    0 0 0 0 mu*Hab_Minus(principal_strains(1),principal_strains(3)) 0; ...
    0 0 0 0 0 mu*Hab_Minus(principal_strains(1),principal_strains(2))];
D_Plus = (A_eps)' * Dn_Plus * A_eps;
D      = A_eps' * (gd*Dn_Plus + Dn_Minus) * A_eps;
D2d = zeros(3,3);
D_Plus_2d = zeros(3,3);
switch nDim
    case 3        
        %转化成张量/矩阵的形式时，需注意：对于2D情况，StressVector（3）代表\sigma_{12}，
        %与应变向量和张量之间的转换关系不一样，应变向量的StrainVector（3）代表2*epsilon_{12}
        Stress_plus = D_Plus * StrainVector;
        sigma_plus = [Stress_plus(1) Stress_plus(6) Stress_plus(5); ...
            Stress_plus(6) Stress_plus(2) Stress_plus(4); ...
            Stress_plus(5) Stress_plus(4) Stress_plus(3)];
    case 2
        D2d = D([1,2,6], [1,2,6]);
        D_Plus_2d = D_Plus([1,2,6], [1,2,6]);
        StressVector = D2d * StrainVector;
        %转化成张量/矩阵的形式时，需注意：对于2D情况，StressVector（3）代表\sigma_{12}，
        %与应变向量和张量之间的转换关系不一样，应变向量的StrainVector（3）代表2*epsilon_{12}
        sigma_plus = D_Plus_2d * StrainVector;
        % sigma_plus = [Stress_plus(1) Stress_plus(3) 0; ...
        %     Stress_plus(3) Stress_plus(2) 0; ...
        %      0 0 0];
         % sigma_plus = [sigma_plus(1,1); sigma_plus(2,2); sigma_plus(1,2)];
end
psi_plus_T = lambda * angle_bracket_Plus(trace(StrainTensor)) * (-2*alpha)...
    + mu*(2 * angle_bracket_Plus(principal_strains(1)) * (-alpha)...
    + 2 * angle_bracket_Plus(principal_strains(2)) * (-alpha)...
    + 2 * angle_bracket_Plus(principal_strains(3)) * (-alpha));%改成sum形式，7.5e-6改成alpha
Psi_negative = lambda/2 * angle_bracket_Minus(trace(StrainTensor))^2 + ...
    mu * sum(angle_bracket_Minus(principal_strains).^2);
k = eps;
Psi = ((1-d)^2 + k) * Psi_plus + Psi_negative;
end
function [s] =Hab_Plus(a,b)
angle_bracket_Plus = @(a) (a+abs(a))/2;
if a==b 
    s=heaviside(a);
else
    s=(angle_bracket_Plus(a) - angle_bracket_Plus(b))/(a-b);
end
end
function [s] =Hab_Minus(a,b)
angle_bracket_Minus = @(a) (a-abs(a))/2;
if a==b 
    s=heaviside(-a);
else
    s=(angle_bracket_Minus(a) - angle_bracket_Minus(b))/(a-b);
end
end
%%
% function [result]=heaviside(x)
% if x>0
%     result=1;
% elseif x==0
%     result=0.5;
% else
%     result=0;
% end



