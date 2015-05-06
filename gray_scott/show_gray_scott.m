%% one proces
clear all
load pro_0.txt -ascii

contourf(pro_0);
title('256^2 points on 1 processes')

%% four processes
clear all
load pro_0.txt -ascii
load pro_1.txt -ascii
load pro_2.txt -ascii
load pro_3.txt -ascii

A=[pro_0 pro_1 ;...
   pro_2 pro_3];

contourf(A)
title('256^2 points on 4 processes')
%% sixteen processes
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
title('256^2 points on 16 processes')
%% thirty-six processes
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
load pro_16.txt -ascii
load pro_17.txt -ascii
load pro_18.txt -ascii
load pro_19.txt -ascii
load pro_20.txt -ascii
load pro_21.txt -ascii
load pro_22.txt -ascii
load pro_23.txt -ascii
load pro_24.txt -ascii
load pro_25.txt -ascii
load pro_26.txt -ascii
load pro_27.txt -ascii
load pro_28.txt -ascii
load pro_29.txt -ascii
load pro_30.txt -ascii
load pro_31.txt -ascii
load pro_32.txt -ascii
load pro_33.txt -ascii
load pro_34.txt -ascii
load pro_35.txt -ascii

A=[pro_0  pro_1  pro_2  pro_3  pro_4  pro_5;...
   pro_6  pro_7  pro_8  pro_9  pro_10 pro_11;...
   pro_12 pro_13 pro_14 pro_15 pro_16 pro_17;...
   pro_18 pro_19 pro_20 pro_21 pro_22 pro_23;...
   pro_24 pro_25 pro_26 pro_27 pro_28 pro_29;...
   pro_30 pro_31 pro_32 pro_33 pro_34 pro_35];

contourf(A)
title('396^2 points on 36 processes')
%% sixty-four processes
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
load pro_16.txt -ascii
load pro_17.txt -ascii
load pro_18.txt -ascii
load pro_19.txt -ascii
load pro_20.txt -ascii
load pro_21.txt -ascii
load pro_22.txt -ascii
load pro_23.txt -ascii
load pro_24.txt -ascii
load pro_25.txt -ascii
load pro_26.txt -ascii
load pro_27.txt -ascii
load pro_28.txt -ascii
load pro_29.txt -ascii
load pro_30.txt -ascii
load pro_31.txt -ascii
load pro_32.txt -ascii
load pro_33.txt -ascii
load pro_34.txt -ascii
load pro_35.txt -ascii
load pro_36.txt -ascii
load pro_37.txt -ascii
load pro_38.txt -ascii
load pro_39.txt -ascii
load pro_40.txt -ascii
load pro_41.txt -ascii
load pro_42.txt -ascii
load pro_43.txt -ascii
load pro_44.txt -ascii
load pro_45.txt -ascii
load pro_46.txt -ascii
load pro_47.txt -ascii
load pro_48.txt -ascii
load pro_49.txt -ascii
load pro_50.txt -ascii
load pro_51.txt -ascii
load pro_52.txt -ascii
load pro_53.txt -ascii
load pro_54.txt -ascii
load pro_55.txt -ascii
load pro_56.txt -ascii
load pro_57.txt -ascii
load pro_58.txt -ascii
load pro_59.txt -ascii
load pro_60.txt -ascii
load pro_61.txt -ascii
load pro_62.txt -ascii
load pro_63.txt -ascii

A=[pro_0  pro_1  pro_2  pro_3  pro_4  pro_5  pro_6  pro_7;...
   pro_8  pro_9  pro_10 pro_11 pro_12 pro_13 pro_14 pro_15;...
   pro_16 pro_17 pro_18 pro_19 pro_20 pro_21 pro_22 pro_23;...
   pro_24 pro_25 pro_26 pro_27 pro_28 pro_29 pro_30 pro_31;...
   pro_32 pro_33 pro_34 pro_35 pro_36 pro_37 pro_38 pro_39;...
   pro_40 pro_41 pro_42 pro_43 pro_44 pro_45 pro_46 pro_47;...
   pro_48 pro_49 pro_50 pro_51 pro_52 pro_53 pro_54 pro_55;...
   pro_56 pro_57 pro_58 pro_59 pro_60 pro_61 pro_62 pro_63];

contourf(A)
title('256^2 points on 64 processes')
%% speedup
T1 = 342;
T4 = 248;
T16 = 59.7;
T36 = 17;

P = [1 4 16 36];

Sp = [T1/T1 T1/T4 T1/T16 T1/T36];

plot(P,Sp, '*')
title('Parallel speedup for 400x400 gray-scott system')
xlabel('Number of processes')
ylabel('Speedup')
