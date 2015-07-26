clear all
load pro_0.txt -ascii
load pro_1.txt -ascii
load pro_2.txt -ascii
load pro_3.txt -ascii

A=[pro_0 pro_1 ;...
   pro_2 pro_3];

contourf(A)
title('300^2 points on 4 processes')