function x = cholsolve(R, b)
%CHOLSOLVE solves the system R'Rx = b using backward- and forward
%substitution
x = mldivide(R, mldivide(R', b));
end