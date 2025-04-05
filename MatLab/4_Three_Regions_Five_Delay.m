%---------------------------------------------------------------------------
% 							Three Regions Five Delay
%---------------------------------------------------------------------------
%Alli3Solver_Region.m
%This is a script file that calls two function files, Alli3Delay_Region,
%and Alli3Delay_hist_region. The purpose of this file is to numerically solve
%the 1 delay alligator model using Matlab"s built in dde solver, Runge
%Kutta Method, dde23.
%Comments updated for AMSSI 2007 on July 27,2007
%%% 2. VARIABLE INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = 0; %starting time value
tF = 200; %the final time for which we will compute a solution
tspan = [t0 tF]; %vector of initial and final time values
b2=.844; %b is the birth rate and is set at a desired value
%%% 1. PARAMETER VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau_1 = 1*b2; %set our delay value tau_1
tau_2 = 8*b2; %set our delay value tau_2
tau_3 = 11*b2; %set our delay value tau_2
tau_4 = 12*b2; %set our delay value tau_2
tau_5 = 31*b2; %set our delay value tau_2
lags = [tau_1 tau_2 tau_3 tau_4 tau_5]; % "lags" is the array that holds the time delay values
%%% 3. NUMERICAL SOLUTION OF ODES %%%%%%%%%%%%%%%%%%%%%%%
%options = ddeset("RelTol",1e-4,"AbsTol",1e-7,"InitialY",[3 5 2 2 4 0 0 0 0 0 0 0 2 2 1 0 0 0 0 0 0 0 0 0]);
AlliBirthDivsol = dde23(@Alli6Delay_Region, lags, @Alli6Delay_hist_region, tspan);% options); %this line calls the
%delay differential equation solver dde23. "Alli6Delay_Region" is the name
%of the file storing out dde formulas, lags is the array of delays that we
%need to use, "Alli6Delay_his_region" is the function file that maintains
%the system history and "tspan" feeds in our initial and final times for
%computation. The solutions to this system will be stored in the array
%"AlliBirthDivsol".
%%% 4. PLOTTING OF SOLUTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%first, separate the output array into its various components
t = AlliBirthDivsol.x; %assigns the time outputs to a vector "t"
Af1 = AlliBirthDivsol.y(1,:); %assigns females between ages 0-1 in R1 outputs
%array to a vector named after the variable
%it represents, "Af1".
Af2 = AlliBirthDivsol.y(2,:); %assigns females between ages 0-1 in R2 outputs
%array to a vector named after the variable
%it represents, "Af2".
Am2 = AlliBirthDivsol.y(3,:); %assigns males between ages 0-1 in R2 outputs array
%to a vector named after the variable it
%represents, "Am2".
Am3 = AlliBirthDivsol.y(4,:); %assigns males between ages 0-1 in R3 outputs array
%to a vector named after the variable it
%represents, "Am3".
Bf1 = AlliBirthDivsol.y(5,:); %assigns females between ages 1-8 in R1 outputs
%array to a vector named after the variable
%it represents, "Bf1".
Bf2 = AlliBirthDivsol.y(6,:); %assigns females between ages 1-8 in R2 outputs
%array to a vector named after the variable
%it represents, "Bf2".
Bm2 = AlliBirthDivsol.y(7,:); %assigns males between ages 1-8 in R2 outputs
%array to a vector named after the variable
%it represents, "Bm2".
Bm3 = AlliBirthDivsol.y(8,:); %assigns males between ages 1-8 in R3 outputs
%array to a vector named after the variable
%it represents, "Bm3".
Cf1 = AlliBirthDivsol.y(9,:); %assigns females between ages 8-11 in R1 outputs
%array to a vector named after the variable
%it represents, "Cf1".
Cf2 = AlliBirthDivsol.y(10,:); %assigns females between ages 8-11 in R2 outputs
%array to a vector named after the variable
%it represents, "Cf2".
Cm2 = AlliBirthDivsol.y(11,:); %assigns males between ages 8-11 in R2 outputs
%array to a vector named after the variable
%it represents, "Cm2".
Cm3 = AlliBirthDivsol.y(12,:); %assigns males between ages 8-11 in R3 outputs
%array to a vector named after the variable
%it represents, "Cm3".
Df1 = AlliBirthDivsol.y(13,:); %assigns females between the ages 11-12 in R1
%outputs array to a vector named after the
%variable it represents, "Df1".
Df2 = AlliBirthDivsol.y(13,:); %assigns females between the ages 11-12 in R2
%outputs array to a vector named after the
%variable it represents, "Df2".
Dm2 = AlliBirthDivsol.y(13,:); %assigns males between the ages 11-12 in R2
%outputs array to a vector named after the
%variable it represents, "Dm2".
Dm3 = AlliBirthDivsol.y(13,:); %assigns males between the ages 11-12 in R3
%outputs array to a vector named after the
%variable it represents, "Dm3".
Ef1 = AlliBirthDivsol.y(17,:); %assigns females between the ages 12-31 in R1
%outputs array to a vector named after the
%variable it represents, "Ef1".
Ef2 = AlliBirthDivsol.y(18,:); %assigns females between the ages 12-31 in R2
%outputs array to a vector named after the
%variable it represents, "Ef2".
Em2 = AlliBirthDivsol.y(19,:); %assigns males between the ages 12-31 in R2
%outputs array to a vector named after the
%variable it represents, "Em2".
Em3 = AlliBirthDivsol.y(20,:); %assigns males between the ages 12-31 in R3
%outputs array to a vector named after the
%variable it represents, "Em3".
Ff1 = AlliBirthDivsol.y(21,:); %assigns females between the ages 31+ in R1
%outputs array to a vector named after the
%variable it represents, "Ff1".
Ff2 = AlliBirthDivsol.y(22,:); %assigns females between the ages 31+ in R2
%outputs array to a vector named after the
%variable it represents, "Ff2".
Fm2 = AlliBirthDivsol.y(23,:); %assigns males between the ages 31+ in R2
%outputs array to a vector named after the
%variable it represents, "Fm2".
Fm3 = AlliBirthDivsol.y(24,:); %assigns males between the ages 31+ in R3
%outputs array to a vector named after the
%variable it represents, "Fm3".
u=Af1+Bf1+Cf1+Df1+Ef1+Ff1; %assigns the variable "u" to be equal to the total female
%population in R1
v=Af2+Bf2+Cf2+Df2+Ef2+Ff2; %assigns the variable "v" to be equal to the total female
%population in R2
w=Am2+Bm2+Cm2+Dm2+Em2+Fm2; %assigns the variable "w" to be equal to the total male
%population in R2
x=Am3+Bm3+Cm3+Dm3+Em3+Fm3; %assigns the variable "x" to be equal to the total male
%population in R3
y=u+v; %assigns the variable "y" to be equal to the total female population
z=w+x; %assigns the variable "z" to be equal to the total male population
r=z./(y+z); %assigns r to be equal to the total male population divided
%by the total male population plus the total female population
%will give you the ratio of males to total population
final=size(t);
Af1(final(2));
Bf1(final(2));
Cf1(final(2));
Df1(final(2));
Ef1(final(2));
Ff1(final(2));
Af2(final(2));
Bf2(final(2));
Cf2(final(2));
Df2(final(2));
Ef2(final(2));
Ff2(final(2));
Am2(final(2));
Bm2(final(2));
Cm2(final(2));
Dm2(final(2));
Em2(final(2));
Fm2(final(2));
Am3(final(2));
Bm3(final(2));
Cm3(final(2));
Dm3(final(2));
Em3(final(2));
Fm3(final(2));
figure % opens a new figure window
hold on % this keeps all the following plots on the same graph
%plot3(t, x, y, "r-", "LineWidth",3) % L will be red and solid of weight 3
plot(t, Af1, "g-", "LineWidth",3) %plot the females between ages 0-1 in R1 a green dashed line
plot(t, Af2, "r-", "LineWidth",3) %plot the females between ages 0-1 in R2 a red dashed line
plot(t, Am2, "y-", "LineWidth",3) %plot the males between ages 0-1 in R2 a yellow dashed line
plot(t, Am3, "b-", "LineWidth",3) %plot the males between ages 0-1 in R3 a blue dashed line
plot(t, Bf1, "g.", "LineWidth",3)%plot the females between ages 1-8 in R1 a green dotted line
plot(t, Bf2, "r.", "LineWidth",3)%plot the females between ages 1-8 in R2 a red dotted line
plot(t, Bm2, "y.", "LineWidth",3)%plot the males between ages 1-8 in R2 a yellow dotted line
plot(t, Bm3, "b.", "LineWidth",3)%plot the males between ages 1-8 in R3 a blue dotted line
plot(t, Cf1, "g-.", "LineWidth",3)%plot the females between ages 8-11 in R1 a green dashed dotted line
plot(t, Cf2, "r-.", "LineWidth",3)%plot the females between ages 8-11 in R2 a red dashed dotted line
plot(t, Cm2, "y-.", "LineWidth",3)%plot the males between ages 8-11 in R2 a yellow dashed dotted line
plot(t, Cm3, "b-.", "LineWidth",3)%plot the males between ages 8-11 in R3 a blue dashed dotted line
plot(t, Df1, "g", "LineWidth",3)%plot the females between ages 11-12 in R1 a green solid line
plot(t, Df1, "r", "LineWidth",3)%plot the females between ages 11-12 in R2 a red solid line
plot(t, Df1, "y", "LineWidth",3)%plot the males between ages 11-12 in R2 a yellow solid line
plot(t, Df1, "b", "LineWidth",3)%plot the males between ages 11-12 in R3 a blue solid line
xlabel("time") % creates the label "time" on the x-axis
ylabel("alligators") %creates the label "stuff" on the y axis
title("Population of ages 0-12") %creates the title
figure %opens a new figure window
hold on %keeps all the following plots on the same graph
plot(t,Ef1,"g-","Linewidth",3); %plot the females between ages 12-31 in R1 a green dashed line
plot(t,Ef2,"r-", "Linewidth",3); %plot the females between ages 12-31 in R2 a red dashed line
plot(t,Em2,"y-","Linewidth",3); %plot the males between ages 12-31 in R2 a yellow dashed line
plot(t,Em3,"b-","Linewidth",3); %plot the males between ages 12-31 in R3 a blue dashed line
plot(t,Ff1,"g","Linewidth",3); %plot the females between ages 31+ in R1 a green solid line
plot(t,Ff2,"r", "Linewidth",3); %plot the females between ages 31+ in R2 a red solid line
plot(t,Fm2,"y","Linewidth",3); %plot the males between ages 31+ in R2 a yellow solid line
plot(t,Fm3,"b","Linewidth",3); %plot the males between ages 31+ in R3 a blue solid line
xlabel("time") % creates the label "time" on the x-axis
ylabel("alligators") %creates the label "stuff" on the y axis
title("Population of ages 12+") %creates the title
figure %open a new figure window
hold on %this keeps all the following plots on the same graph
plot(t, r,"m","Linewidth",3) %plot the ratio of males to the total population a red line
title("the love graph") %creates the title
%total_f1_pde=3070.45;
%total_f2_pde=466.385;
%total_m2_pde=252.79;
%total_m3_pde=233.705;
%bargraphvector1=[(u(final(2))/(u(final(2))+v(final(2))+w(final(2))+x(final(2))))...
% (v(final(2))/(u(final(2))+v(final(2))+w(final(2))+x(final(2))))...
% (w(final(2))/(u(final(2))+v(final(2))+w(final(2))+x(final(2))))...
% (x(final(2))/(u(final(2))+v(final(2))+w(final(2))+x(final(2))))]
%bargraphvector2=[total_f1_pde/(total_f1_pde+total_f2_pde+total_m2_pde+total_m3_pde)...
% total_f2_pde/(total_f1_pde+total_f2_pde+total_m2_pde+total_m3_pde)...
% total_m2_pde/(total_f1_pde+total_f2_pde+total_m2_pde+total_m3_pde)...
% total_m3_pde/(total_f1_pde+total_f2_pde+total_m2_pde+total_m3_pde)]
%figure
%bar([bargraphvector1; bargraphvector2],"group")
%colormap spring
%legend("F1","F2","M2","M3")
% creates the label "time" on the x-axis
%ylabel("percentage of total population") %creates the label "stuff" on the y axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a function code that stores all the delay differential equations
%formulas in another file to be called by the script file
function dydt=Alli6Delay_Region(t,y,Z) %function name and variables, time
% variables and lags
%Current size(Z) =[3,2]
%col 1 of Z: solutions with time delay
%tau_1. col 2 of Z represents solutions
%with time delay tau_2.
%row 1 of Z: delay value for x
%row 2 of Z: delay value for y
%row 3 of Z: delay value for y
global d1 d2 d3 b1 b2 s1 s8 s11 s12 s31 k1 k2 k3 c1 c2 r1 a1 a2 a3
%declare variables as global to allow the dde23 solver to access their
%values, but we only have to change them in one place. The global command
%is to be used before you have defined the variables in the code.
d1=.9; %d1 is the death rate between the ages 0-8 and is set at a desired value
d2=.151; %d2 is the death rate between the ages 8-11 and is set at a desired value
d_new=.001;%d_new is the death rate between the ages 11-31 and is set a desired value
%we gave the age group 11-31 a very low death rate since biological data
%states that it was 0.
d3=.139; %d3 is the death rate between the ages 31+ and is set at a desired value
b1=.286; %b1 is the birth rate for the teens and is set at a desired value
b2=.844; %b2 is the birth rate for the adults and is set at a desired value
k1=.79*1.5; %k1 is the carrying capacity of R1 and is set at a desired value
k2=.15*1.5; %k2 is the carrying capacity of R2 and is set at a desired value
k3=.6*1.5; %k3 is the carrying capacity of R3 and is set at a desired value
s1=exp(-d1);%s1 is the survival probability to age 1 and is set at a desired value
s8=exp(-d2.*7);%s8 is the survival probability from age 1 to age 8 and is set at a desired value
s11 =exp(-d2.*3); %s11 is the survival probability from age 8 to age 11 and is set at a desired value
s12 =exp(-d_new); %s12 is the survival probability from age 11 to age 12 and is set at a desired value
s31 =exp(-d_new*19);%s31 is the survival probability from age 12 to age 31+ and is set at a desired value
c1=k2/k1; %c1 is the carrying capacity of R2 divided by the carrying capacity of R1
c2=k3/k1; %c2 is the carrying capacity of R3 divided by the carrying capacity of R1
r1=b1/b2; %r1 is the birth rate for the teens divided by the birth rate for adults
a1=d1./b2; %a1 is the death rate between the ages 0-8 divided by the birth rate for the adults
a2=d2./b2; %a2 is the death rate between the ages 8-11 divided by the birth rate for the adults
a_new=d_new/b2; %a_new is the death rate between the ages 11-31 divided by the birth rate for the adults
a3=d3./b2; %a3 is the death rate between the ages 31+ divided by the birth rate for the adults
Lag1=Z(:,1); %a column vector whose rows represent the solution at the time delay 1
Lag2=Z(:,2); %a column vector whose rows represent the solution at the time delay 2
Lag3=Z(:,3); %a column vector whose rows represent the solution at the time delay 3
Lag4=Z(:,4); %a column vector whose rows represent the solution at the time delay 4
Lag5=Z(:,5); %a column vector whose rows represent the solution at the time delay 5
Q1=(y(9)+y(13)+y(17)+y(21))/(1+y(9)+y(13)+y(17)+y(21));
%this is the fraction of f1 who can"t nest in R1 now
Q2=(c1/(c1+y(9)+y(13)+y(17)+y(21)+y(10)+y(14)+y(18)+y(22)));
%this is the fraction of f1s and f2s who can nest in R2 now
Q3=(Lag1(9)+Lag1(13)+Lag1(17)+Lag1(21))/(1+Lag1(9)+Lag1(13)+Lag1(17)+Lag1(21));
%this is the fraction of f1 who couldn"t nest in R1 1 year ago
Q4=(c1/(c1+Lag1(9)+Lag1(13)+Lag1(17)+Lag1(21)+Lag1(10)+Lag1(14)+Lag1(18)+Lag1(22)));
%this is the fraction of f1 and f2 who could nest in R2 1 year ago
Q5=(c2/(c2+y(9)+y(13)+y(17)+y(21)+y(10)+y(14)+y(18)+y(22)));
%this is the fraction of f1s and f2s who can nest in R3 now
Q6=(c2/(c2+Lag1(9)+Lag1(13)+Lag1(17)+Lag1(21)+Lag1(10)+Lag1(14)+Lag1(18)+Lag1(22)));
%this is the fraction of f1s and f2s who could nest in R2 1 year ago
R3=(Lag2(9)+Lag2(13)+Lag2(17)+Lag2(21))/(1+Lag2(9)+Lag2(13)+Lag2(17)+Lag2(21));
%this is the fraction of f1 who couldn"t nest in R1 8 years ago
R4=(c1/(c1+Lag2(9)+Lag2(13)+Lag2(17)+Lag2(21)+Lag2(10)+Lag2(14)+Lag2(18)+Lag2(22)));
%this is the fraction of f1 and f2 who could nest in R2 8 years ago
R6=(c2/(c2+Lag2(9)+Lag2(13)+Lag2(17)+Lag2(21)+Lag2(10)+Lag2(14)+Lag2(18)+Lag2(22)));
%this is the fraction of f1s and f2s who could nest in R2 8 years ago
S3=(Lag3(9)+Lag3(13)+Lag3(17)+Lag3(21))/(1+Lag3(9)+Lag3(13)+Lag3(17)+Lag3(21));
%this is the fraction of f1 who couldn"t nest in R1 11 years ago
S4=(c1/(c1+Lag3(9)+Lag3(13)+Lag3(17)+Lag3(21)+Lag3(10)+Lag3(14)+Lag3(18)+Lag3(22)));
%this is the fraction of f1 and f2 who could nest in R2 11 years ago
S6=(c2/(c2+Lag3(9)+Lag3(13)+Lag3(17)+Lag3(21)+Lag3(10)+Lag3(14)+Lag3(18)+Lag3(22)));
%this is the fraction of f1s and f2s who could nest in R3 11 years ago
T3=(Lag4(9)+Lag4(13)+Lag4(17)+Lag4(21))/(1+Lag4(9)+Lag4(13)+Lag4(17)+Lag4(21));
%this is the fraction of f1 who couldn"t nest in R1 12 years ago
T4=(c1/(c1+Lag4(9)+Lag4(13)+Lag4(17)+Lag4(21)+Lag4(10)+Lag4(14)+Lag4(18)+Lag4(22)));
%this is the fraction of f1 and f2 who could nest in R2 12 years ago
T6=(c2/(c2+Lag4(9)+Lag4(13)+Lag4(17)+Lag4(21)+Lag4(10)+Lag4(14)+Lag4(18)+Lag4(22)));
%this is the fraction of f1s and f2s who could nest in R2 12 years ago
U3=(Lag5(9)+Lag5(13)+Lag5(17)+Lag5(21))/(1+Lag5(9)+Lag5(13)+Lag5(17)+Lag5(21));
%this is the fraction of f1 who couldn"t nest in R1 31 years ago
U4=(c1/(c1+Lag5(9)+Lag5(13)+Lag5(17)+Lag5(21)+Lag5(10)+Lag5(14)+Lag5(18)+Lag5(22)));
%this is the fraction of f1 and f2 who could nest in R2 31 years ago
U6=(c2/(c2+Lag5(9)+Lag5(13)+Lag5(17)+Lag5(21)+Lag5(10)+Lag5(14)+Lag5(18)+Lag5(22)));
%this is the fraction of f1s and f2s who could nest in R2 31 years ago
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 0-1 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 0-1:y(1), R2 females between 0-1:y(2), R2 males between 0-1:y(3), R3 males between 0-1:y(4)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(1)=r1*((y(9)+y(13))*(1-Q1))-s1.*r1*(Lag1(9)+Lag1(13))*(1-Q3)...
-a1*y(1)+...
(y(17)+y(21))*(1-Q1)...
-s1.*(Lag1(17)+Lag1(21))*(1-Q3);
%delay differential equation for the rate of change of the females between
%ages 0-1 R1 where Cf1(t-tau_1) is represented by Lag1(9), Df1(t-tau_1) is
%represented by Lag1(13), Ef1(t-tau_1) is represented by Lag1(17),and
%Ff1(t-tau_1) is represented by Lag1(21).
dydt(2)=(r1/2)*((y(10)+y(14)+(y(9)+y(13))*Q1)*Q2...
-s1.*(Lag1(10)+Lag1(14)+(Lag1(9)+Lag1(13))*Q3)*Q4)...
-a1*y(2)...
+(1/2)*((y(18)+y(22)+((y(17)+y(21))*Q1)*Q2)...
-s1.*(Lag1(18)+Lag1(22)+(Lag1(17)+Lag1(21))*Q3)*Q4);
%delay differential equation for the rate of change of the females between
%ages 0-1 R2 where Cf2(t-tau_1) is represented by Lag1(10), Df2(t-tau_1) is
%represented by Lag1(14), Cf1(t-tau_1) is represented by Lag1(9),
%Df1(t-tau_1) is represented by Lag1(13), Ef2(t-tau_1) is represented by
%Lag1(18), Ff2(t-tau_1) is represented by Lag1(22), Ef1(t-tau_1) is
%represented by Lag1(17),Ff1(t-tau_1) is represented by Lag1(22).
dydt(3)=(r1/2)*((y(10)+y(14)+(y(9)+y(13))*Q1)*Q2...
-s1.*(Lag1(10)+Lag1(14)+(Lag1(9)+Lag1(13))*Q3)*Q4)...
-a1*y(3)...
+(1/2)*((y(18)+y(22)+((y(17)+y(21))*Q1)*Q2)...
-s1.*(Lag1(18)+Lag1(22)+(Lag1(17)+Lag1(21))*Q3)*Q4);
%delay differential equation for the rate of change of the males between
%ages 0-1 R3 where Cf2(t-tau_1) is represented by Lag1(10), Df2(t-tau_1) is
%represented by Lag1(14), Cf1(t-tau_1) is represented by Lag1(9),
%Df1(t-tau_1) is represented by Lag1(13), Ef2(t-tau_1) is represented by
%Lag1(18), Ff2(t-tau_1) is represented by Lag1(22), Ef1(t-tau_1) is
%represented by Lag1(17),Ff1(t-tau_1) is represented by Lag1(22).
dydt(4)=r1*(Q5*((y(9)+y(13))*Q1+y(10)+y(14))*(1-Q2)...
-s1.*Q6*((Lag1(9)+Lag1(13))*Q3+Lag1(10)+Lag1(14))*(1-Q4))...
-a1*y(4)...
+Q5*((y(17)+y(21))*Q1+y(18)+y(22))*(1-Q2)...
-s1.*Q6*((Lag1(17)+Lag1(21))*Q3+Lag1(18)+Lag1(22))*(1-Q4);
%delay differential equation for the rate of change of the males between
%ages 0-1 R3 where Cf1(t-tau_1) is represented by Lag1(9), Df1(t-tau_1) is
%represented by Lag1(13), Cf2(t-tau_1) is represented by Lag1(10),
%Df2(t-tau_1) is represented by Lag1(14), Ef1(t-tau_1) is represented by
%Lag1(17), Ff1(t-tau_1) is represented by Lag1(21), Ef2(t-tau_1) is
%represented by Lag1(18),Ff1(t-tau_1) is represented by Lag1(22).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 1-8 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 1-8:y(5), R2 females between 1-8:y(6), R2 males between 1-8:y(7), R3 males between 1-8:y(8)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(5)=s1.*r1*(Lag1(9)+Lag1(13))*(1-Q3)+...
s1.*(Lag1(17)+Lag1(21))*(1-Q3)...
-a2*y(5)-...
s1.*s8.*r1*(Lag2(9)+Lag2(13))*(1-R3)-...
s1.*s8.*(Lag2(17)+Lag2(21))*(1-R3);
%delay differential equation for the rate of change of the females between
%ages 1-8 R1 where Cf1(t-tau_1) is represented by Lag1(9), Df1(t-tau_1) is
%represented by Lag1(13), Ef1(t-tau_1) is represented by Lag1(17),
%Ff1(t-tau_1) is represented by Lag1(21), Cf1(t-tau_2) is represented by
%Lag2(9), Df1(t-tau_2) is represented by Lag2(13), Ef1(t-tau_2) is
%represented by Lag2(17), and Ff1(t-tau_2) is represented by Lag2(21).
dydt(6)=(r1/2)*s1.*(Lag1(10)+Lag1(14)+(Lag1(9)+Lag1(13))*Q3)*Q4+...
(1/2)*s1.*(Lag1(18)+Lag1(22)+(Lag1(17)+Lag1(21))*Q3)*Q4...
-a2*y(6)-...
(r1/2)*s1.*s8.*(Lag2(10)+Lag2(14)+(Lag2(9)+Lag2(13))*R3)*R4-...
(1/2)*s1.*s8.*(Lag2(18)+Lag2(22)+(Lag2(17)+Lag2(21))*R3)*R4;
%delay differential equation for the rate of change of the females between
%ages 1-8 R2 where Cf2(t-tau_1) is represented by Lag1(10), Df2(t-tau_1) is
%represented by Lag1(14), Cf1(t-tau_1) is represented by Lag1(9),
%Df1(t-tau_1) is represented by Lag1(13), Ef2(t-tau_1) is represented by
%Lag1(18), Ff2(t-tau_1) is represented by Lag1(22), Ef1(t-tau_1) is
%represented by Lag1(17), Ff1(t-tau_1) is represented by Lag1(21),
%Cf2(t-tau_2) is represented by Lag2(10), Df2(t-tau_2) is
%represented by Lag2(14), Cf1(t-tau_2) is represented by Lag2(9),
%Df1(t-tau_2) is represented by Lag2(13), Ef2(t-tau_2) is represented by
%Lag2(18), Ff2(t-tau_2) is represented by Lag2(22), Ef1(t-tau_2) is
%represented by Lag2(17), and Ff1(t-tau_2) is represented by Lag2(21).
dydt(7)=(r1/2)*s1.*(Lag1(10)+Lag1(14)+(Lag1(9)+Lag1(13))*Q3)*Q4+...
(1/2)*s1.*(Lag1(18)+Lag1(22)+(Lag1(17)+Lag1(21))*Q3)*Q4...
-a2*y(7)-...
(r1/2)*s1.*s8.*(Lag2(10)+Lag2(14)+(Lag2(9)+Lag2(13))*R3)*R4-...
(1/2)*s1.*s8.*(Lag2(18)+Lag2(22)+(Lag2(17)+Lag2(21))*R3)*R4;
%delay differential equation for the rate of change of the males between
%ages 1-8 R2 where Cf2(t-tau_1) is represented by Lag1(10), Df2(t-tau_1) is
%represented by Lag1(14), Cf1(t-tau_1) is represented by Lag1(9),
%Df1(t-tau_1) is represented by Lag1(13), Ef2(t-tau_1) is represented by
%Lag1(18), Ff2(t-tau_1) is represented by Lag1(22), Ef1(t-tau_1) is
%represented by Lag1(17), Ff1(t-tau_1) is represented by Lag1(21),
%Cf2(t-tau_2) is represented by Lag2(10), Df2(t-tau_2) is
%represented by Lag2(14), Cf1(t-tau_2) is represented by Lag2(9),
%Df1(t-tau_2) is represented by Lag2(13), Ef2(t-tau_2) is represented by
%Lag2(18), Ff2(t-tau_2) is represented by Lag2(22), Ef1(t-tau_2) is
%represented by Lag2(17), and Ff1(t-tau_2) is represented by Lag2(21).
dydt(8)=r1*s1.*Q6*((Lag1(9)+Lag1(13))*Q3+Lag1(10)+Lag1(14))*(1-Q4)+...
s1.*Q6*((Lag1(17)+Lag1(21))*Q3+Lag1(18)+Lag1(22))*(1-Q4)-...
a2*y(8)-...
r1*s1.*s8.*R6*((Lag2(9)+Lag2(13))*R3+Lag2(10)+Lag2(14))*(1-R4)-...
s1.*s8.*R6*((Lag2(17)+Lag2(21))*R3+Lag2(18)+Lag2(22))*(1-R4);
%delay differential equation for the rate of change of the males between
%ages 1-8 R3 where Cf1(t-tau_1) is represented by Lag1(9), Df1(t-tau_1) is
%represented by Lag1(13), Cf2(t-tau_1) is represented by Lag1(10),
%Df2(t-tau_1) is represented by Lag1(14), Ef1(t-tau_1) is represented by
%Lag1(17), Ff1(t-tau_1) is represented by Lag1(21), Ef2(t-tau_1) is
%represented by Lag1(18), Ff2(t-tau_1) is represented by Lag1(22),
%Cf1(t-tau_2) is represented by Lag2(9), Df1(t-tau_2) is
%represented by Lag2(13), Cf2(t-tau_2) is represented by Lag2(10),
%Df2(t-tau_2) is represented by Lag2(14), Ef1(t-tau_2) is represented by
%Lag2(17), Ff1(t-tau_2) is represented by Lag2(21), Ef2(t-tau_2) is
%represented by Lag2(18), and Ff2(t-tau_2) is represented by Lag2(22).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 8-11 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 8-11:y(9), R2 females between 8-11:y(10), R2 males between 8-11:y(11), R3 males between 8-11:y(12)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(9)= s1.*s8.*r1*(Lag2(9)+Lag2(13))*(1-R3)+...
s1.*s8.*(Lag2(17)+Lag2(21))*(1-R3)...
-a2*y(9)-...
s1.*s8.*s11.*r1*(Lag3(9)+Lag3(13))*(1-S3)-...
s1.*s8.*s11.*(Lag3(17)+Lag3(21))*(1-S3);
%delay differential equation for the rate of change of the females between
%ages 8-11 R1 where Cf1(t-tau_2) is represented by Lag2(9), Df1(t-tau_2) is
%represented by Lag2(13), Ef1(t-tau_2) is represented by Lag2(17),
%Ff1(t-tau_2) is represented by Lag2(21), Cf1(t-tau_3) is represented by
%Lag3(9), Df1(t-tau_3) is represented by Lag3(13), Ef1(t-tau_3) is
%represented by Lag3(17), and Ff1(t-tau_3) is represented by Lag3(21).
dydt(10)=(r1/2)*s1.*s8.*(Lag2(10)+Lag2(14)+(Lag2(9)+Lag2(13))*R3)*R4+...
(1/2)*s1.*s8.*(Lag2(18)+Lag2(22)+(Lag2(17)+Lag2(21))*R3)*R4...
-a2*y(10)...
-(r1/2)*s1.*s8.*s11.*(Lag3(10)+Lag3(14)+(Lag3(9)+Lag3(13))*S3)*S4-...
(1/2)*s1.*s8.*s11.*(Lag3(18)+Lag3(22)+(Lag3(17)+Lag3(21))*S3)*S4;
%delay differential equation for the rate of change of the females between
%ages 8-11 R2 where Cf2(t-tau_2) is represented by Lag2(10), Df2(t-tau_2) is
%represented by Lag2(14), Cf1(t-tau_2) is represented by Lag2(9),
%Df1(t-tau_2) is represented by Lag2(13), Ef2(t-tau_2) is represented by
%Lag2(18), Ff2(t-tau_2) is represented by Lag2(22), Ef1(t-tau_2) is
%represented by Lag2(17), Ff1(t-tau_2) is represented by Lag2(21),
%Cf2(t-tau_3) is represented by Lag3(10), Df2(t-tau_3) is
%represented by Lag3(14), Cf1(t-tau_3) is represented by Lag3(9),
%Df1(t-tau_3) is represented by Lag3(13), Ef2(t-tau_3) is represented by
%Lag3(18), Ff2(t-tau_3) is represented by Lag3(22), Ef1(t-tau_3) is
%represented by Lag3(17), and Ff1(t-tau_3) is represented by Lag3(21).
dydt(11)=(r1/2)*s1.*s8.*(Lag2(10)+Lag2(14)+(Lag2(9)+Lag2(13))*R3)*R4+...
(1/2)*s1.*s8.*(Lag2(18)+Lag2(22)+(Lag2(17)+Lag2(21))*R3)*R4...
-a2*y(11)...
-(r1/2)*s1.*s8.*s11.*(Lag3(10)+Lag3(14)+(Lag3(9)+Lag3(13))*S3)*S4-...
(1/2)*s1.*s8.*s11.*(Lag3(18)+Lag3(22)+(Lag3(17)+Lag3(21))*S3)*S4;
%delay differential equation for the rate of change of the males between
%ages 8-11 R2 where Cf2(t-tau_2) is represented by Lag2(10), Df2(t-tau_2) is
%represented by Lag2(14), Cf1(t-tau_2) is represented by Lag2(9),
%Df1(t-tau_2) is represented by Lag2(13), Ef2(t-tau_2) is represented by
%Lag2(18), Ff2(t-tau_2) is represented by Lag2(22), Ef1(t-tau_2) is
%represented by Lag2(17), Ff1(t-tau_2) is represented by Lag2(21),
%Cf2(t-tau_3) is represented by Lag3(10), Df2(t-tau_3) is
%represented by Lag3(14), Cf1(t-tau_3) is represented by Lag3(9),
%Df1(t-tau_3) is represented by Lag3(13), Ef2(t-tau_3) is represented by
%Lag3(18), Ff2(t-tau_3) is represented by Lag3(22), Ef1(t-tau_3) is
%represented by Lag3(17), and Ff1(t-tau_3) is represented by Lag3(21).
dydt(12)=r1*s1.*s8.*R6*((Lag2(9)+Lag2(13))*R3+Lag2(10)+Lag2(14))*(1-R4)+...
s1.*s8.*R6*((Lag2(17)+Lag2(21))*R3+Lag2(18)+Lag2(22))*(1-R4)...
-a2*y(12)...
-r1*s1.*s8.*s11.*S6*((Lag3(9)+Lag3(13))*S3+Lag3(10)+Lag3(14))*(1-S4)...
-s1.*s8.*s11.*S6*((Lag3(17)+Lag3(21))*S3+Lag3(18)+Lag3(22))*(1-S4);
%delay differential equation for the rate of change of the males between
%ages 8-11 R3 where Cf1(t-tau_2) is represented by Lag2(9), Df1(t-tau_2) is
%represented by Lag2(13), Cf2(t-tau_2) is represented by Lag2(10),
%Df2(t-tau_2) is represented by Lag2(14), Ef1(t-tau_2) is represented by
%Lag2(17), Ff1(t-tau_2) is represented by Lag2(21), Ef2(t-tau_2) is
%represented by Lag2(18), Ff2(t-tau_2) is represented by Lag2(22),
%Cf1(t-tau_3) is represented by Lag3(9), Df1(t-tau_3) is
%represented by Lag3(13), Cf2(t-tau_3) is represented by Lag3(10),
%Df2(t-tau_3) is represented by Lag3(14), Ef1(t-tau_3) is represented by
%Lag3(17), Ff1(t-tau_3) is represented by Lag3(21), Ef2(t-tau_3) is
%represented by Lag3(18), and Ff2(t-tau_3) is represented by Lag3(22).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 11-12 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 11-12:y(13), R2 females between 11-12:y(14), R2 males between 11-12:y(15), R3 males between 11-12:y(16)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(13)=s1.*s8.*s11.*r1*(Lag3(9)+Lag3(13))*(1-S3)+...
s1.*s8.*s11.*(Lag3(17)+Lag3(21))*(1-S3)...
-a_new*y(13)...
-s1.*s8.*s11.*s12*r1*(Lag4(9)+Lag4(13))*(1-T3)-...
s1.*s8.*s11.*s12*(Lag4(17)+Lag4(21))*(1-T3);
%delay differential equation for the rate of change of the females between
%ages 11-12 R1 where Cf1(t-tau_3) is represented by Lag3(9), Df1(t-tau_3) is
%represented by Lag3(13), Ef1(t-tau_3) is represented by Lag3(17),
%Ff1(t-tau_3) is represented by Lag3(21), Cf1(t-tau_4) is represented by
%Lag4(9), Df1(t-tau_4) is represented by Lag4(13), Ef1(t-tau_4) is
%represented by Lag4(17), and Ff1(t-tau_4) is represented by Lag4(21).
dydt(14)=(r1/2)*s1.*s8.*s11.*(Lag3(10)+Lag3(14)+(Lag3(9)+Lag3(13))*S3)*S4+...
(1/2)*s1.*s8.*s11.*(Lag3(18)+Lag3(22)+(Lag3(17)+Lag3(21))*S3)*S4...
-a_new*y(14)...
-(r1/2)*s1.*s8.*s11.*s12*(Lag4(10)+Lag4(14)+(Lag4(9)+Lag4(13))*T3)*T4-...
(1/2)*s1.*s8.*s11.*s12*(Lag4(18)+Lag4(22)+(Lag4(17)+Lag4(21))*T3)*T4;
%delay differential equation for the rate of change of the females between
%ages 11-12 R2 where Cf2(t-tau_3) is represented by Lag3(10), Df2(t-tau_3) is
%represented by Lag3(14), Cf1(t-tau_3) is represented by Lag3(9),
%Df1(t-tau_3) is represented by Lag3(13), Ef2(t-tau_3) is represented by
%Lag3(18), Ff2(t-tau_3) is represented by Lag3(22), Ef1(t-tau_3) is
%represented by Lag3(17), Ff1(t-tau_3) is represented by Lag3(21),
%Cf2(t-tau_4) is represented by Lag4(10), Df2(t-tau_4) is
%represented by Lag4(14), Cf1(t-tau_4) is represented by Lag4(9),
%Df1(t-tau_4) is represented by Lag4(13), Ef2(t-tau_4) is represented by
%Lag4(18), Ff2(t-tau_4) is represented by Lag4(22), Ef1(t-tau_4) is
%represented by Lag4(17), and Ff1(t-tau_4) is represented by Lag4(21).
dydt(15)=(r1/2)*s1.*s8.*s11.*(Lag3(10)+Lag3(14)+(Lag3(9)+Lag3(13))*S3)*S4+...
(1/2)*s1.*s8.*s11.*(Lag3(18)+Lag3(22)+(Lag3(17)+Lag3(21))*S3)*S4...
-a_new*y(15)...
-(r1/2)*s1.*s8.*s11.*s12*(Lag4(10)+Lag4(14)+(Lag4(9)+Lag4(13))*T3)*T4-...
(1/2)*s1.*s8.*s11.*s12*(Lag4(18)+Lag4(22)+(Lag4(17)+Lag4(21))*T3)*T4;
%delay differential equation for the rate of change of the males between
%ages 11-12 R2 where Cf2(t-tau_3) is represented by Lag3(10), Df2(t-tau_3) is
%represented by Lag3(14), Cf1(t-tau_3) is represented by Lag3(9),
%Df1(t-tau_3) is represented by Lag3(13), Ef2(t-tau_3) is represented by
%Lag3(18), Ff2(t-tau_3) is represented by Lag3(22), Ef1(t-tau_3) is
%represented by Lag3(17), Ff1(t-tau_3) is represented by Lag3(21),
%Cf2(t-tau_4) is represented by Lag4(10), Df2(t-tau_4) is
%represented by Lag4(14), Cf1(t-tau_4) is represented by Lag4(9),
%Df1(t-tau_4) is represented by Lag4(13), Ef2(t-tau_4) is represented by
%Lag4(18), Ff2(t-tau_4) is represented by Lag4(22), Ef1(t-tau_4) is
%represented by Lag4(17), and Ff1(t-tau_4) is represented by Lag4(21).
dydt(16)=r1*s1.*s8.*s11.*S6*((Lag3(9)+Lag3(13))*S3+Lag3(10)+Lag3(14))*(1-S4)+...
s1.*s8.*s11.*S6*((Lag3(17)+Lag3(21))*S3+Lag3(18)+Lag3(22))*(1-S4)...
-a_new*y(16)...
-r1*s1.*s8.*s11.*s12*T6*((Lag4(9)+Lag4(13))*T3+Lag4(10)+Lag4(14))*(1-T4)...
-s1.*s8.*s11.*s12*T6*((Lag4(17)+Lag4(21))*T3+Lag4(18)+Lag4(22))*(1-T4);
%delay differential equation for the rate of change of the males between
%ages 11-12 R3 where Cf1(t-tau_3) is represented by Lag3(9), Df1(t-tau_3) is
%represented by Lag3(13), Cf2(t-tau_3) is represented by Lag3(10),
%Df2(t-tau_3) is represented by Lag3(14), Ef1(t-tau_3) is represented by
%Lag3(17), Ff1(t-tau_3) is represented by Lag3(21), Ef2(t-tau_3) is
%represented by Lag3(18), Ff2(t-tau_3) is represented by Lag3(22),
%Cf1(t-tau_4) is represented by Lag4(9), Df1(t-tau_4) is
%represented by Lag4(13), Cf2(t-tau_4) is represented by Lag4(10),
%Df2(t-tau_4) is represented by Lag4(14), Ef1(t-tau_4) is represented by
%Lag4(17), Ff1(t-tau_4) is represented by Lag4(21), Ef2(t-tau_4) is
%represented by Lag4(18), and Ff2(t-tau_4) is represented by Lag4(22).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 12-31 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 12-31:y(17), R2 females between 12-31:y(18), R2 males between 12-31:y(19), R3 males between 12-31:y(20)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(17)=s1.*s8.*s11.*s12*r1*(Lag4(9)+Lag4(13))*(1-T3)+...
s1.*s8.*s11.*s12*(Lag4(17)+Lag4(21))*(1-T3)...
-a_new*y(17)...
-s1.*s8.*s11.*s12*s31*r1*(Lag5(9)+Lag5(13))*(1-U3)-...
s1.*s8.*s11.*s12*s31*(Lag5(17)+Lag5(21))*(1-U3);
%delay differential equation for the rate of change of the females between
%ages 12-31 R1 where Cf1(t-tau_4) is represented by Lag4(9), Df1(t-tau_4) is
%represented by Lag4(13), Ef1(t-tau_4) is represented by Lag4(17),
%Ff1(t-tau_4) is represented by Lag4(21), Cf1(t-tau_5) is represented by
%Lag5(9), Df1(t-tau_5) is represented by Lag5(13), Ef1(t-tau_5) is
%represented by Lag5(17), and Ff1(t-tau_5) is represented by Lag5(21).
dydt(18)=(r1/2)*s1.*s8.*s11.*s12*(Lag4(10)+Lag4(14)+(Lag4(9)+Lag4(13))*T3)*T4+...
(1/2)*s1.*s8.*s11.*s12*(Lag4(18)+Lag4(22)+(Lag4(17)+Lag4(21))*T3)*T4-...
a_new*y(18)...
-(r1/2)*s1.*s8.*s11.*s12*s31*(Lag5(10)+Lag5(14)+(Lag5(9)+Lag5(13))*U3)*U4-...
(1/2)*s1.*s8.*s11.*s12*s31*(Lag5(18)+Lag5(22)+(Lag5(17)+Lag5(21))*U3)*U4;
%delay differential equation for the rate of change of the females between
%ages 12-31 R2 where Cf2(t-tau_4) is represented by Lag4(10), Df2(t-tau_4) is
%represented by Lag4(14), Cf1(t-tau_4) is represented by Lag4(9),
%Df1(t-tau_4) is represented by Lag4(13), Ef2(t-tau_4) is represented by
%Lag4(18), Ff2(t-tau_4) is represented by Lag4(22), Ef1(t-tau_4) is
%represented by Lag4(17), Ff1(t-tau_4) is represented by Lag4(21),
%Cf2(t-tau_5) is represented by Lag5(10), Df2(t-tau_5) is
%represented by Lag5(14), Cf1(t-tau_5) is represented by Lag5(9),
%Df1(t-tau_5) is represented by Lag5(13), Ef2(t-tau_5) is represented by
%Lag5(18), Ff2(t-tau_5) is represented by Lag5(22), Ef1(t-tau_5) is
%represented by Lag5(17), and Ff1(t-tau_5) is represented by Lag5(21).
dydt(19)=(r1/2)*s1.*s8.*s11.*s12*(Lag4(10)+Lag4(14)+(Lag4(9)+Lag4(13))*T3)*T4+...
(1/2)*s1.*s8.*s11.*s12*(Lag4(18)+Lag4(22)+(Lag4(17)+Lag4(21))*T3)*T4-...
a_new*y(19)...
-(r1/2)*s1.*s8.*s11.*s12*s31*(Lag5(10)+Lag5(14)+(Lag5(9)+Lag5(13))*U3)*U4-...
(1/2)*s1.*s8.*s11.*s12*s31*(Lag5(18)+Lag5(22)+(Lag5(17)+Lag5(21))*U3)*U4;
%delay differential equation for the rate of change of the males between
%ages 12-31 R2 where Cf2(t-tau_4) is represented by Lag4(10), Df2(t-tau_4) is
%represented by Lag4(14), Cf1(t-tau_4) is represented by Lag4(9),
%Df1(t-tau_4) is represented by Lag4(13), Ef2(t-tau_4) is represented by
%Lag4(18), Ff2(t-tau_4) is represented by Lag4(22), Ef1(t-tau_4) is
%represented by Lag4(17), Ff1(t-tau_4) is represented by Lag4(21),
%Cf2(t-tau_5) is represented by Lag5(10), Df2(t-tau_5) is
%represented by Lag5(14), Cf1(t-tau_5) is represented by Lag5(9),
%Df1(t-tau_5) is represented by Lag5(13), Ef2(t-tau_5) is represented by
%Lag5(18), Ff2(t-tau_5) is represented by Lag5(22), Ef1(t-tau_5) is
%represented by Lag5(17), and Ff1(t-tau_5) is represented by Lag5(21).
dydt(20)=r1*s1.*s8.*s11.*s12*T6*((Lag4(9)+Lag4(13))*T3+Lag4(10)+Lag4(14))*(1-T4)...
+s1.*s8.*s11.*s12*T6*((Lag4(17)+Lag4(21))*T3+Lag4(18)+Lag4(22))*(1-T4)...
-a_new*y(20)...
-r1*s1.*s8.*s11.*s12*s31*U6*((Lag5(9)+Lag5(13))*U3+Lag5(10)+Lag5(14))*(1-U4)...
-s1.*s8.*s11.*s12*s31*U6*((Lag5(17)+Lag5(21))*U3+Lag5(18)+Lag5(22))*(1-U4);
%delay differential equation for the rate of change of the males between
%ages 12-31 R3 where Cf1(t-tau_4) is represented by Lag4(9), Df1(t-tau_4) is
%represented by Lag4(13), Cf2(t-tau_4) is represented by Lag4(10),
%Df2(t-tau_4) is represented by Lag4(14), Ef1(t-tau_4) is represented by
%Lag4(17), Ff1(t-tau_4) is represented by Lag4(21), Ef2(t-tau_4) is
%represented by Lag4(18), Ff2(t-tau_4) is represented by Lag4(22),
%Cf1(t-tau_5) is represented by Lag5(9), Df1(t-tau_5) is
%represented by Lag5(13), Cf2(t-tau_5) is represented by Lag5(10),
%Df2(t-tau_5) is represented by Lag5(14), Ef1(t-tau_5) is represented by
%Lag5(17), Ff1(t-tau_5) is represented by Lag5(21), Ef2(t-tau_5) is
%represented by Lag5(18), and Ff2(t-tau_5) is represented by Lag5(22).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 31+ Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 31+:y(21), R2 females between 31+:y(22), R2 males between 31+:y(23), R3 males between 31+:y(24)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(21)=s1.*s8.*s11.*s12*s31*r1*(Lag5(9)+Lag5(13))*(1/(1+Lag5(9)+Lag5(13)+Lag5(17)+Lag5(21)))+...
s1.*s8.*s11.*s12*s31*(Lag5(17)+Lag5(21))*(1/(1+Lag5(9)+Lag5(13)+Lag5(17)+Lag5(21)))...
-a3*y(21);
%delay differential equation for the rate of change of the females between
%ages 31+ R1 where Cf1(t-tau_5) is represented by Lag5(9), Df1(t-tau_5) is
%represented by Lag5(13), Ef1(t-tau_5) is represented by Lag5(17), and
%Ff1(t-tau_5) is represented by Lag5(21).
dydt(22)=(r1/2)*s1.*s8.*s11.*s12*s31*(Lag5(10)+Lag5(14)+(Lag5(9)+Lag5(13))*U3)*U4+...
(1/2)*s1.*s8.*s11.*s12*s31*(Lag5(18)+Lag5(22)+(Lag5(17)+Lag5(21))*U3)*U4...
-a3*y(22);
%delay differential equation for the rate of change of the females between
%ages 31+ R2 where Cf2(t-tau_5) is represented by Lag5(10), Df2(t-tau_5)
%is represented by Lag5(14), Cf1(t-tau_5) is represented by Lag5(9),
%Df1(t-tau_5) is represented by Lag5(13), Ef2(t-tau_5) is represented by
%Lag5(18), Ff2(t-tau_5) is represented by Lag5(22), Ef1(t-tau_5) is
%represented by Lag5(17), and Ff1(t-tau_5) is represented by Lag5(21).
dydt(23)=(r1/2)*s1.*s8.*s11.*s12*s31*(Lag5(10)+Lag5(14)+(Lag5(9)+Lag5(13))*U3)*U4+...
(1/2)*s1.*s8.*s11.*s12*s31*(Lag5(18)+Lag5(22)+(Lag5(17)+Lag5(21))*U3)*U4...
-a3*y(23);
%delay differential equation for the rate of change of the males between
%ages 31+ R2 where Cf2(t-tau_5) is represented by Lag5(10), Df2(t-tau_5)
%is represented by Lag5(14), Cf1(t-tau_5) is represented by Lag5(9),
%Df1(t-tau_5) is represented by Lag5(13), Ef2(t-tau_5) is represented by
%Lag5(18), Ff2(t-tau_5) is represented by Lag5(22), Ef1(t-tau_5) is
%represented by Lag5(17), and Ff1(t-tau_5) is represented by Lag5(21).
dydt(24)=r1*s1.*s8.*s11.*s12*s31*U6*((Lag5(9)+Lag5(13))*U3+Lag5(10)+Lag5(14))*(1-U4)...
+s1.*s8.*s11.*s12*s31*U6*((Lag5(17)+Lag5(21))*U3+Lag5(18)+Lag5(22))*(1-U4)...
-a3*y(24);
%delay differential equation for the rate of change of the males between
%ages 31+ R3 where Cf1(t-tau_5) is represented by Lag5(9), Df1(t-tau_5)
%is represented by Lag5(13), Cf2(t-tau_5) is represented by Lag5(10),
%Df2(t-tau_5) is represented by Lag5(14), Ef1(t-tau_5) is represented by
%Lag5(17), Ff1(t-tau_5) is represented by Lag5(21), Ef2(t-tau_5) is
%represented by Lag5(18), and Ff2(t-tau_5) is represented by Lag5(22).
dydt=dydt'; %the solver must return a column vector and this sets Matlab
%up to return a row vector, so we take the transpose to keep it
%all the same.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%DDE23 requires histories for all variables even if they do not have a delay in
%the equation. Due to evidence the death die out after a certain amount of time we choose
%constants for our histories.
%This is a function file that maintains the system histories.
function s = Alli6Delay_hist_region(t) %function name and input
%for dde23, the output needs to be a column
%vector
%In this case we are assigning constant values to the variable histories.
hist_1 = 5; %this is the history of the females between ages 0-1 population from R1
hist_2 = 5; %this is the history of the females between ages 0-1 population from R2
hist_3 = 5; %this is the history of the males between ages 0-1 population from R2
hist_4 = 5; %this is the history of the males between ages 0-1 population from R3
hist_5 = 5; %this is the history of the females between ages 1-8 population from R1
hist_6 = 5; %this is the history of the females between ages 1-8 population from R2
hist_7 = 5; %this is the history of the males between ages 1-8 population from R2
hist_8 = 5; %this is the history of the males between ages 1-8 population from R3
hist_9 = 5; %this is the history of the females between ages 8-11 population from R1
hist_10 = 3.5; %this is the history of the females between ages 8-11 population from R2
hist_11 = 2; %this is the history of the males between ages 8-11 population from R2
hist_12 =2; %this is the history of the males between ages 8-11 population from R3
hist_13 = 5; %this is the history of the females between ages 11-12 population from R1
hist_14 = 3.5; %this is the history of the females between ages 11-12 population from R2
hist_15 = 2; %this is the history of the males between ages 11-12 population from R2
hist_16 = 2; %this is the history of the males between ages 11-12 population from R3
hist_17 = 2.5; %this is the history of the females between ages 12-31 population from R1
hist_18 = 1.5; %this is the history of the females between ages 12-31 population from R2
hist_19 = 1; %this is the history of the males between ages 12-31 population from R2
hist_20 = 1; %this is the history of the males between ages 12-31 population from R3
hist_21 = 2.5; %this is the history of the females between ages 31+ population from R1
hist_22 = 1.5; %this is the history of the females between ages 31+ population from R2
hist_23 = 1; %this is the history of the males between ages 31+ population from R2
hist_24 = 1; %this is the history of the males between ages 31+ population from R3
s = [hist_1; hist_2; hist_3; hist_4; hist_5; hist_6; hist_7;hist_8;hist_9;hist_10;hist_11;hist_12;
hist_13;hist_14;hist_15;hist_16;hist_17;hist_18;hist_19;hist_20;hist_21;hist_22;hist_23;hist_24];
%This puts all the histories into a column vector and names it s. This is
%the format that dde23 would like it to be in.
end