uudd=sqrt(max(abs(eigs((Kd_active\Kdu_active)*(Ku_active\Kud_active)))))
ddtt=sqrt(max(abs(eigs((KT_active\KTd_active)*(Kd_active\KdT_active)))))
uutt=sqrt(max(abs(eigs((KT_active\zeros(size(KuT_active))')*(Ku_active\KuT_active)))))
%uutt=sqrt(max(abs(eigs((KT_active\ones(size(KuT_active))')*(Ku_active\ones(size(KuT_active)))))))
