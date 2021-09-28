function vec_out_action=action(Vr,Threshold,V_peak,dt,time_of_action,t_start,Time_vec)

Time_vec1=Time_vec(1):round(mean(Time_vec));
Time_vec2=Time_vec1(end)+1:Time_vec(end);

A1=(V_peak-Threshold)./(Time_vec1(end)-Time_vec1(1));
B1=Threshold-(A1*Time_vec1(1));
V1=A1.*Time_vec1+B1;


A2=(Vr-V_peak)./(Time_vec2(end)-Time_vec2(1));
B2=V_peak-(A2*Time_vec2(1));
V2=A2.*Time_vec2+B2;




vec_out_action=[V1 V2];




