% test QSVT in solving QLSP
% ----------------------------------------------------------------------
%
% Author:           Yulong Dong, dongyl@berkeley.edu
% Version:          1.0
% Last revision:    5/11/2020
%
%  ----------------------------------------------------------------------

clear
clc

% parameters
n_sys_qubit = 5;
herm_perturb = false;

N = 2^n_sys_qubit;
scale_fac = [];
n_sample = 2000;
degodd = zeros(3,1);
precision = zeros(3,1);
dec_rat = 10^(-8/n_sample);

for kappa = 10:10:30
    % inverse general matrix needs odd QSVT
    asset = "odd"+int2str(kappa);
    data = load("phasefactor_xinv.mat",asset);
    assignin("base","phi",data.(asset));
    asset = "scale_fac"+int2str(kappa);
    data = load("phasefactor_xinv.mat",asset);
    assignin("base","scale_fac",data.(asset));
    degodd(round(kappa/10)) = length(phi)-1;

    % initiate a complex matrix with controllable condition number
    [W,~] = qr(rand(N,N)+1j*rand(N,N));
    if herm_perturb
        V = W;
    else
        [V,~] = qr(rand(N,N)+1j*rand(N,N));
    end
    A = W * diag(rand(N,1)*(1-1/kappa)+1/kappa) * V';
    % Block-encode A by SVD
    U = BlockEncode_Easy(A);
    if ~herm_perturb
        % interchange U, U'
        U.mat = U.mat';
    end

    retoddqc =  QSVT_QC(U,phi,true);
    % compute exact matrix after scaling
    Ainv = inv(A) / scale_fac;

    fprintf("QSVT kappa: %d, Odd: %5e\n",kappa, norm(retoddqc - Ainv))
    precision(round(kappa/10)) = norm(retoddqc - Ainv);
    % retnorm, perturbnorm
    result = zeros(n_sample,2);
    perturb_eps = 1e-2;
    for jj = 1:n_sample
        [delta_W,~] = qr(rand(N,N)+1j*rand(N,N));
        if herm_perturb
            delta_V = delta_W;
        else
            [delta_V,~] = qr(rand(N,N)+1j*rand(N,N));
        end
        delta_A = perturb_eps * delta_W * diag(rand(N,1)) * delta_V';
        U_perturb = BlockEncode_Easy(A + delta_A);
        if ~herm_perturb
            U_perturb.mat = U_perturb.mat';
        end
        retoddqc_perturb =  QSVT_QC(U_perturb,phi,true);
        norm_odd = norm(retoddqc_perturb-retoddqc);
        norm_perturb = norm(delta_A);
        result(jj,:) = [norm_odd, norm_perturb];
        perturb_eps = perturb_eps * dec_rat;
    end
    assignin("base","odd"+int2str(kappa),result);
end
if herm_perturb
    save("data_qsvt_herm.mat","odd10","odd20","odd30","degodd","precision","-v7.3");
else
    save("data_qsvt.mat","odd10","odd20","odd30","degodd","precision","-v7.3");
end

