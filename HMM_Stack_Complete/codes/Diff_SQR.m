function [Diff] = Diff_SQR(old_core_param,core_param)

Diff = 0;

for p = 1:length(core_param)
    Diff = Diff + (old_core_param(p).shift - core_param(p).shift)*(old_core_param(p).shift - core_param(p).shift);
    Diff = Diff + (old_core_param(p).R - core_param(p).R)*(old_core_param(p).R - core_param(p).R);
    Diff = Diff + (old_core_param(p).alpha - core_param(p).alpha)*(old_core_param(p).alpha - core_param(p).alpha);
    Diff = Diff + (old_core_param(p).beta - core_param(p).beta)*(old_core_param(p).beta - core_param(p).beta);
    Diff = Diff + (old_core_param(p).phi - core_param(p).phi)*(old_core_param(p).phi - core_param(p).phi);
    Diff = Diff + (old_core_param(p).psi - core_param(p).psi)*(old_core_param(p).psi - core_param(p).psi);
    Diff = Diff + (old_core_param(p).eta - core_param(p).eta)*(old_core_param(p).eta - core_param(p).eta);
    Diff = Diff + (old_core_param(p).epsilon - core_param(p).epsilon)*(old_core_param(p).epsilon - core_param(p).epsilon);
end

Diff = sqrt(Diff/length(core_param));

end

