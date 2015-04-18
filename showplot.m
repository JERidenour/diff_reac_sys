load u_vals.txt -ascii

u = u_vals;
n = sqrt(length(u));
contourf(reshape(u,n,n))