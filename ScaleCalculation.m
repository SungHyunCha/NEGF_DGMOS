function result = ScaleCalculation(A,b,VT,ni)

[X,Y] = size(A);
C = sparse(1:3:X-2,1:3:Y-2,VT,X,Y)...
    +sparse(2:3:X-1,2:3:Y-1,ni,X,Y)...
    +sparse(3:3:X,3:3:Y,ni,X,Y);
AC = sparse(A*C);
R = sum(abs(AC),2).^-1;
R = spdiags(R,0,X,Y);
result = C*((R*AC)\(R*b));

end