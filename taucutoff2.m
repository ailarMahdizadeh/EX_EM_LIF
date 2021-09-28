clc
close all
clear;
n=1;%n number of vesicles
g_ca=12*10^(-12);

Vr = -50;Threshold = 10;V_peak=75;
time_of_action=2;
RefPeriod =2; %ms
%dimension equalization
dimension_k_pos=1e+9;dimension_k_neg=1e+3;
P=1.6;RTF=(26.7)*10^(-3);Caex=10*10^(-3);A=0.1*10^(9);
kplus1=3.75*10^(-3)*dimension_k_pos;kminus1=4*10^(-4)*dimension_k_neg;
kplus2=2.5*10^(-3)*dimension_k_pos;kminus2=10^(-3)*dimension_k_neg;
kplus3=5*10^(-4)*dimension_k_pos;kminus3=0.1*dimension_k_neg;
kplus4=7.5*10^(-3)*dimension_k_pos;kminus4=10*dimension_k_neg;
kplus=[kplus1 kplus2 kplus3 kplus4];
kminus=[kminus1 kminus2 kminus3 kminus4];
end_time=50; %ms
dt = 0.5;
time = 0:dt:end_time;
num_Ref_period=RefPeriod./dt;
start_of_freq=0.1; %Hz
%% frequency initialization
end_of_freq=500; %Hz
step_of_freq=end_of_freq/100; %Hz
freq=start_of_freq:step_of_freq:end_of_freq;
Prel=zeros(numel(time),numel(freq));
spike_number=zeros(numel(freq),numel(time));

R=100;G=1/R;
%% tau change
%  Markers = {'+','o','*','x','v','d','^','s','>','<'};
% set(0,'defaultaxeslinestyleorder',{'-',':','o','*','x','v','d','^','s','>','<'}) %or whatever you want
set(groot,'defaultaxeslinestyleorder','remove')
ax.FontSize=17;
tau_check=[1  2  3 4 5]; %ms
V_m=zeros(numel(freq),numel(tau_check));
for index_tau=1:numel(tau_check)
    tau_LIF=tau_check(index_tau);
    C=tau_LIF*G;
for j=1:numel(freq)
j
%%
T_pulse=1e3/freq(j); %ms
pulsewidth = T_pulse/2.2; %ms
d=[0:T_pulse:end_time]';
xx = @rectpuls;
%% Excitory
% I0=pulstran(time,d,xx,pulsewidth);
 I0=(sin(2*pi*freq(j).*time.*1e-3));

%%  main
V = zeros(size(time));
V(1) = Vr;
t = 2;
while t<size(time,2)
   
      if V(t-1)>=Threshold
         t;
          spike_number(j,t-1)=(sum(spike_number(j,1:t-2))+1)/((t-1)*dt);
          Time_vec= t:t+(time_of_action./dt);
          V(Time_vec) =action(Vr,Threshold,V_peak,dt,time_of_action,t,Time_vec); %vec_out_action=action(Vr,Threshold,V_peak,dt,time_of_action,t_start,Time_vec)
    
     V(Time_vec(end)+1:Time_vec(end)+num_Ref_period)=Vr;
      t=Time_vec(end)+num_Ref_period;
      else
          V(t) = V(t-1) + dt*(Vr-V(t-1) + R.*I0(t))/tau_LIF;
          t=t+1; 
      end 
     
  
end
V(numel(time)+1:end)=[];
V=V.*1e-3;
V_m(j,index_tau)=max(V);

end
legendInfo{index_tau} = ['\tau = ' num2str(tau_check(index_tau)) ' [ms]']; 


%% plot

hold on;
% plot(freq,V_m(:,index_tau)*1e3,'LineWidth',2,Markers{index_tau})
% plot(freq,V_m(:,index_tau)*1e3,'LineWidth',2,strcat('-',Markers{index_tau}))
 %set(gca,'LineStyle',Markers(index_tau)) 
 if index_tau==1
plot(freq,V_m(:,index_tau)*1e3,'LineWidth',2)
 elseif  index_tau==2
plot(freq,V_m(:,index_tau)*1e3,'LineWidth',2,'LineStyle',':')
 elseif index_tau==3
   plot(freq,V_m(:,index_tau)*1e3,'LineWidth',2,'LineStyle','-.')   
   elseif index_tau==4
   plot(freq,V_m(:,index_tau)*1e3,'-o','LineWidth',2)   
 else
  plot(freq,V_m(:,index_tau)*1e3,'LineWidth',2,'LineStyle','--') 
 end
 
end
grid on
xlabel('Excitation Frequency [Hz]')
ylabel('Maximum membrane voltage [mV]')
legend(legendInfo)
