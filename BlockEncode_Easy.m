function [mat_block_encode] = BlockEncode_Easy(mat)
% given a square matrix mat whose 2-norm is bounded by 1
% output a (1,1,0)-block-encoding of it by SVD
% signal state is |0\rangle
% ----------------------------------------------------------------------
%
% Author:           Yulong Dong, dongyl@berkeley.edu
% Version:          1.0
% Last revision:    5/11/2020
%
%  ----------------------------------------------------------------------

[row, col] = size(mat);
if row ~= col
    error("BlockEncode_Easy: matrix must be square")
end
[W, S, V] = svd(mat);

if max(S) > 1
    error("BlockEncode_Easy: matrix 2-norm must be bounded by 1");
end

Ssqrt = sqrt(eye(row) - S^2);
ret = [mat, W*Ssqrt; Ssqrt*V', -S];
mat_block_encode.mat = ret;
mat_block_encode.m = 1;
mat_block_encode.n = round(log2(row));

end

