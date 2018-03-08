function [ matrix ] = vectorToMatrix( vector, nx, ny )
matrix = reshape(vector, [nx, ny]);
end

