function [ vector ] = matrixToVector( matrix, nx, ny )
vector = reshape(matrix, [nx*ny, 1]);
end
