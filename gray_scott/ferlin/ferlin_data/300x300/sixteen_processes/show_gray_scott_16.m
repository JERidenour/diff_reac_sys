clear all
load pro_0.txt -ascii
load pro_1.txt -ascii
load pro_2.txt -ascii
load pro_3.txt -ascii
load pro_4.txt -ascii
load pro_5.txt -ascii
load pro_6.txt -ascii
load pro_7.txt -ascii
load pro_8.txt -ascii
load pro_9.txt -ascii
load pro_10.txt -ascii
load pro_11.txt -ascii
load pro_12.txt -ascii
load pro_13.txt -ascii
load pro_14.txt -ascii
load pro_15.txt -ascii

A=[pro_0  pro_1  pro_2  pro_3;...
   pro_4  pro_5  pro_6  pro_7;...
   pro_8  pro_9  pro_10 pro_11;...
   pro_12 pro_13 pro_14 pro_15];

contourf(A)
title('300^2 points on 16 processes')