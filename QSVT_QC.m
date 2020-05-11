function [ret] = QSVT_QC(blockencode, phi, real_part)
% QSVT circuit for generic block encoding
% blockencode - unitary matrix block-encodes a matrix with bounded norm
% phi - phase factor in SU(2)-QSVT
% real_part - bool
% ret - f(A) with global phase
% ----------------------------------------------------------------------
%
% Author:           Yulong Dong, dongyl@berkeley.edu
% Version:          1.0
% Last revision:    5/11/2020
%
%  ----------------------------------------------------------------------

m = blockencode.m;
n = blockencode.n;
% turn phi to \varphi
phase = phi + pi/2;
phase(1) = phase(1) - pi/4;
phase(end) = phase(end) - pi/4;

sigmaz = diag([1,-1]);
sigmax = [0,1;1,0];
Hadamard = [1,1;1,-1] / sqrt(2);
Gstate = sparse([1],[1],[1],2^m,1);
Gproj = Gstate*Gstate';
CNOT = kron(eye(2), speye(2^m) - Gproj) + kron(sigmax, Gproj);
CNOT = kron(CNOT, speye(2^n));
UA = kron(eye(2), blockencode.mat);
Zgate = kron(sigmaz, speye(2^(n+m)));
Hgate = kron(Hadamard, speye(2^(n+m)));

if real_part
    ret = Hgate;
    if mod(length(phase)-1,2) == 1
        ret = Zgate * ret;
    end
else
    ret = speye(2^(n+m+1));
end
dag = false;
for jj = length(phase):-1:2
    rotmat = [exp(-1i*phase(jj)),0;0,exp(1i*phase(jj))];
    rotmat = kron(rotmat, speye(2^(n+m)));
    if dag
        ret = UA' * CNOT * rotmat * CNOT * ret;
    else
        ret = UA * CNOT * rotmat * CNOT * ret;
    end
    dag = ~dag;
end
rotmat = [exp(-1i*phase(1)),0;0,exp(1i*phase(1))];
rotmat = kron(rotmat, speye(2^(n+m)));
ret = CNOT * rotmat * CNOT * ret;
if real_part
    ret = Hgate * ret;
end

% add global phase, Eq.(12)
ret = (-1i)^(length(phase)-1)*ret;

ret = ret(1:2^n, 1:2^n);

end

