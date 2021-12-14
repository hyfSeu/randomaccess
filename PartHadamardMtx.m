function [Phi] = PartHadamardMtx(M, N)
    L_t = max(M,N);
    L_t1 = (12 - mod(L_t,12)) + L_t;
    L_t2 = (20 - mod(L_t,20)) + L_t;
    L_t3 = 2^ceil(log2(L_t));
    L = min([L_t1,L_t2,L_t3]);
    Phi = [];
    Phi_t = hadamard(L);
    RowIndex = randperm(L);
    Phi_t_r = Phi_t(RowIndex(1:M),:);
    ColIndex = randperm(L);
    Phi = Phi_t_r(:,ColIndex(1:N));
end