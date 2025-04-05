%---------------------------------------------------------------------------
% 							Three Regions Three Delay
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
b2=.844; %%b is the birth rate and is set at a desired value
%%% 1. PARAMETER VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tau_1 = 1*b2; %set our delay value tau_1
tau_2 = 8*b2; %set our delay value tau_2
tau_3 = 11*b2; %set our delay value tau_2
%tau_4 = 12*b2; %set our delay value tau_2
tau_5 = 31*b2; %set our delay value tau_2
lags = [tau_2 tau_3 tau_5]; %"lags" is the array that holds the time delay values
%%% 3. NUMERICAL SOLUTION OF ODES %%%%%%%%%%%%%%%%%%%%%%%
AlliBirthDivsol = dde23(@Alli3Delay_Region, lags, @Alli3Delay_hist_region, tspan);
%delay differential equation solver dde23. "Alli3Delay_Region" is the name
%of the file storing out dde formulas, lags is the array of delays that we
%need o use, "Alli3Delay_his_region" is the function file that maintains
%the system history and "tspan" feeds in our initial and final times for
%computation. The solutions to this system will be stored in the array
%"AlliBirthDivsol".
%%% 4. PLOTTING OF SOLUTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first, separate the output array into its various components
t = AlliBirthDivsol.x; %assigns the time outputs to a vector "t"
Af1 = AlliBirthDivsol.y(1,:); %assigns females between ages 0-8 in R1 outputs
%array to a vector named after the variable
%it represents, "Af1".
Af2 = AlliBirthDivsol.y(2,:); %assigns females between ages 0-8 in R2 outputs
%array to a vector named after the variable
%it represents, "Af2".
Am2 = AlliBirthDivsol.y(3,:); %assigns males between ages 0-8 in R2 outputs array
%to a vector named after the variable it
%represents, "Am2".
Am3 = AlliBirthDivsol.y(4,:); %assigns males between ages 0-8 in R3 outputs array
%to a vector named after the variable it
%represents, "Am3".
Bf1 = AlliBirthDivsol.y(5,:); %assigns females between ages 8-11 in R1 outputs
%array to a vector named after the variable
%it represents, "Bf1".
Bf2 = AlliBirthDivsol.y(6,:); %assigns females between ages 8-11 in R2 outputs
%array to a vector named after the variable
%it represents, "Bf2".
Bm2 = AlliBirthDivsol.y(7,:); %assigns males between ages 8-11 in R2 outputs
%array to a vector named after the variable
%it represents, "Bm2".
Bm3 = AlliBirthDivsol.y(8,:); %assigns males between ages 8-11 in R3 outputs
%array to a vector named after the variable
%it represents, "Bm3".
Cf1 = AlliBirthDivsol.y(9,:); %assigns females between ages 11-31 in R1 outputs
%array to a vector named after the variable
%it represents, "Cf1".
Cf2 = AlliBirthDivsol.y(10,:); %assigns females between ages 11-31 in R2 outputs
%array to a vector named after the variable
%it represents, "Cf2".
Cm2 = AlliBirthDivsol.y(11,:); %assigns males between ages 11-31 in R2 outputs
%array to a vector named after the variable
%it represents, "Cm2".
Cm3 = AlliBirthDivsol.y(12,:); %assigns males between ages 11-31 in R3 outputs
%array to a vector named after the variable
%it represents, "Cm3".
Df1 = AlliBirthDivsol.y(13,:); %assigns females between the ages 31+ in R1
%outputs array to a vector named after the
%variable it represents, "Df1".
Df2 = AlliBirthDivsol.y(13,:); %assigns females between the ages 31+ in R2
%outputs array to a vector named after the
%variable it represents, "Df2".
Dm2 = AlliBirthDivsol.y(13,:); %assigns males between the ages 31+ in R2
%outputs array to a vector named after the
%variable it represents, "Dm2".
Dm3 = AlliBirthDivsol.y(13,:); %assigns males between the ages 31+ in R3
%outputs array to a vector named after the
%variable it represents, "Df1".
u=Af1+Bf1+Cf1+Df1; %assigns the variable "u" to be equal to the total female
%population in R1
v=Af2+Bf2+Cf2+Df2; %assigns the variable "v" to be equal to the total female
%population in R2
w=Am2+Bm2+Cm2+Dm2; %assigns the variable "w" to be equal to the total male
%population in R2
q=Am3+Bm3+Cm3+Dm3; %assigns the variable "q" to be equal to the total male
%population in R3
y=u+v; %assigns the variable "y" to be equal to the total female population
z=w+q; %assigns the variable "z" to be equal to the total male population
r=z./(y+z); %assigns r to be equal to the total male population divided
%by the total male population plus the total female population
%will give you the ratio of males to total population
figure % opens a new figure window
hold on % this keeps all the following plots on the same graph
%plot3(t, x, y, "r-", "LineWidth",3) % L will be red and solid of weight 3
plot(t, Af1, "g-", "LineWidth",3) %plot the females between ages 0-8 in R1 a green dashed line
plot(t, Af2, "r-", "LineWidth",3) %plot the females between ages 0-8 in R2 a red dashed line
plot(t, Am2, "y-", "LineWidth",3) %plot the males between ages 0-8 in R2 a yellow dashed line
plot(t, Am3, "b-", "LineWidth",3) %plot the males between ages 0-8 in R3 a blue dashed line
plot(t, Bf1, "g.", "LineWidth",3)%plot the females between ages 8-11 in R1 a green dotted line
plot(t, Bf2, "r.", "LineWidth",3)%plot the females between ages 8-11 in R2 a red dotted line
plot(t, Bm2, "y.", "LineWidth",3)%plot the males between ages 8-11 in R2 a yellow dotted line
plot(t, Bm3, "b.", "LineWidth",3)%plot the males between ages 8-11 in R3 a blue dotted line
plot(t, Cf1, "g-.", "LineWidth",3)%plot the females between ages 11-31 in R1 a green dashed dotted line
plot(t, Cf2, "r-.", "LineWidth",3)%plot the females between ages 11-31 in R2 a red dashed dotted line
plot(t, Cm2, "y-.", "LineWidth",3)%plot the males between ages 11-31 in R2 a yellow dashed dotted line
plot(t, Cm3, "b-.", "LineWidth",3)%plot the males between ages 11-31 in R3 a blue dashed dotted line
plot(t, Df1, "g", "LineWidth",3)%plot the females between ages 31+ in R1 a green solid line
plot(t, Df1, "r", "LineWidth",3)%plot the females between ages 31+ in R2 a red solid line
plot(t, Df1, "y", "LineWidth",3)%plot the males between ages 31+ in R2 a yellow solid line
plot(t, Df1, "b", "LineWidth",3)%plot the males between ages 31+ in R3 a blue solid line
xlabel("time") % creates the label "time" on the x-axis
ylabel("alligators") %creates the label "stuff" on the y axis
title("F1 (m) and F2 (c) and M2 (y) M3 (b) J1(solid) J2(dotted) A(dashed)") %creates the title
figure %open a new figure window
hold on %this keeps all the following plots on the same graph
plot(t, u, "g", "LineWidth",3) %plot the total female population in R1 a green line
plot(t, v, "r", "LineWidth",3) %plot the total female population in R2 a red line
plot(t, w, "y", "LineWidth",3) %plot the total male population in R2 a yellow line
plot(t, q, "b", "LineWidth",3) %plot the total male population in R3 a blue line
figure %open a new figure window
hold on %this keeps all the following plots on the same graph
plot(t,r,"r","Linewidth",3) %plot the ratio of males to the total population a red line
title("the love graph") %creates the title
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a function code that stores all the delay differential equations
%formulas in another file to be called by the script file
function dydt=Alli3Delay_Region(t,y,Z) %function name and variables, time
% variables and lags
global d1 d2 d3 b1 b2 s8 s11 k1 k2 k3 c1 c2 r1 a1 a2 a3
%declare variables as global to allow the dde23 solver to access their
%values, but we only have to change them in one place. The global command
%is to be used before you have defined the variables in the code.
d1=.245; %d1 is the death rate between the ages 0-8 and is set at a desired value
d2=.151; %d2 is the death rate between the ages 8-11 and is set at a desired value
d_new=.001; %d_new is the death rate between the ages 11-31 and is set a desired value
%we gave the age group 11-31 a very low death rate since biological data
%states that it was 0.
d3=.139; %d3 is the death rate between the ages 31+ and is set at a desired value
b1=.286; %b1 is the birth rate for the teens and is set at a desired value
b2=.844; %bA is the birth rate for the adults and is set at a desired value
k1=.79*1.5; %k1 is the carrying capacity of R1 and is set at a desired value
k2=.15*1.5; %k2 is the carrying capacity of R2 and is set at a desired value
k3=.06*1.5; %k3 is the carrying capacity of R3 and is set at a desired value
s8=exp(-d1*8); %s8 is the survival probability to age 8 and is set at a desired value
s11 =exp(-d2*3);%s11 is the survival probability from age 8 to age 11 and is set at a desired value
s31=exp(-d_new*20); %s31 is the survival probability from age 11 to age 31 and is set at a desired value
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
Q1=(y(5)+y(9)+y(13))/(1+y(5)+y(9)+y(13));
%this is the fraction of f1 who can"t nest in R1 now
Q2=(c1/(c1+y(5)+y(6)+y(9)+y(10)+y(13)+y(14)));
%this is the fraction of f1s and f2s who can nest in R2 now
Q3=(Lag1(5)+Lag1(9)+Lag1(13))/(1+Lag1(5)+Lag1(9)+Lag1(13));
%this is the fraction of f1 who couldn"t nest in R1 8 years ago
Q4=(c1/(c1+Lag1(5)+Lag1(6)+Lag1(9)+Lag1(10)+Lag1(13)+Lag1(14)));
%this is the fraction of f1 and f2 who could nest in R2 8 years ago
Q5=(c2/(c2+y(5)+y(6)+y(9)+y(10)+y(13)+y(14)));
%this is the fraction of f1s and f2s who can nest in R3 now
Q6=(c2/(c2+Lag1(5)+Lag1(6)+Lag1(9)+Lag1(10)+Lag1(13)+Lag1(14)));
%this is the fraction of f1s and f2s who could nest in R2 8 years ago
R3=(Lag2(5)+Lag2(9)+Lag2(13))/(1+Lag2(5)+Lag2(9)+Lag2(13));
%this is the fraction of f1 who couldn"t nest in R1 11 years ago
R4=(c1/(c1+Lag2(5)+Lag2(6)+Lag2(9)+Lag2(10)+Lag2(13)+Lag2(14)));
%this is the fraction of f1 and f2 who could nest in R2 11 years ago
R6=(c2/(c2+Lag2(5)+Lag2(6)+Lag2(9)+Lag2(10)+Lag2(13)+Lag2(14)));
%this is the fraction of f1s and f2s who could nest in R2 11 years ago
S3=(Lag3(5)+Lag3(9)+Lag3(13))/(1+Lag3(5)+Lag3(9)+Lag3(13));
%this is the fraction of f1 who couldn"t nest in R1 11 years ago
S4=(c1/(c1+Lag3(5)+Lag3(6)+Lag3(9)+Lag3(10)+Lag3(13)+Lag3(14)));
%this is the fraction of f1 and f2 who could nest in R2 11 years ago
S6=(c2/(c2+Lag3(5)+Lag3(6)+Lag3(9)+Lag3(10)+Lag3(13)+Lag3(14)));
%this is the fraction of f1s and f2s who could nest in R2 11 years ago
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 0-8 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 0-8:y(1), R2 females between 0-8:y(2), R2 males between 0-8:y(3), R3 males between 0-8:y(4)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(1)=r1*(y(5)*(1-Q1))...
-s8*r1*(Lag1(5))*(1-Q3)...
-a1*y(1)...
+(y(9)+y(13))*(1-Q1)...
-s8*(Lag1(9)+Lag1(13))*(1-Q3);
%delay differential equation for the rate of change of the females between
%ages 0-8 R1 where Bf1(t-tau_2) is represented by Lag1(5), Cf1(t-tau_2) is
%represented by Lag1(9), and Df1(t-tau_2) is represented by Lag1(13).
dydt(2)=(r1/2)*(((y(5)+y(6))*Q1)*Q2...
-s8*(Lag1(5)+Lag1(6))*Q3*Q4)...
-a1*y(2)...
+(1/2)*((y(9)+y(10)+y(13)+y(14))*Q1)*Q2...
-s8*(Lag1(9)+Lag1(10)+Lag1(13)+Lag1(14))*Q3*Q4;
%delay differential equation for the rate of change of the females between
%ages 0-8 in R2 where Bf1(t-tau_2) is represented by Lag1(5), Bf2(t-tau_2)
%is represented by Lag1(6),Cf1(t-tau_2) is represented by Lag1(9),
%Cf2(t-tau_2) is represented by Lag1(10), Df1(t-tau_2) is represnted by
%Lag1(13), and Df1(t-tau_2) is represented by Lag1(14).
dydt(3)=(r1/2)*(((y(5)+y(6))*Q1)*Q2...
-s8*(Lag1(5)+Lag1(6))*Q3*Q4)...
-a1*y(3)...
+(1/2)*((y(9)+y(10)+y(13)+y(14))*Q1)*Q2...
-s8*(Lag1(9)+Lag1(10)+Lag1(13)+Lag1(14))*Q3*Q4;
%delay differential equation for the rate of change of the males between
%ages 0-8 R2 where Bf1(t-tau_2) is represented by Lag1(5), Bf2(t-tau_2) is
%represented by Lag1(6),Cf1(t-tau_2) is represented by Lag1(9),
%Cf2(t-tau_2) is represented by Lag1(10), Df1(t-tau_2) is represnted by
%Lag1(13), and Df1(t-tau_2) is represented by Lag1(14).
dydt(4)=r1*(Q5*(y(5)*Q1+y(6))*(1-Q2)...
-s8*Q6*(Lag1(5)*Q3+Lag1(6))*(1-Q4))...
-a1*y(4)...
+Q5*((y(9)+y(13))*Q1+y(10)+y(14))*(1-Q2)...
-s8*Q6*((Lag1(9)+Lag1(13))*Q3+Lag1(10)+Lag1(14))*(1-Q4);
%delay differential equation for the rate of change of the males between
%ages 0-8 R3 where Bf1(t-tau_2) is represented by Lag1(5), Bf2(t-tau_2) is
%represented by Lag1(6),Cf1(t-tau_2) is represented by Lag1(9),
%Cf2(t-tau_2) is represented by Lag1(10), Df1(t-tau_2) is represnted by
%Lag1(13), and Df1(t-tau_2) is represented by Lag1(14).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 8-11 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 8-11:y(5), R2 females between 8-11:y(6), R2 males between 8-11:y(7), R3 males between 8-11:y(8)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(5)=s8*r1*(Lag1(5))*(1-Q3)...
+s8*(Lag1(9)+Lag1(13))*(1-Q3)...
-a2*y(5)...
-s8*s11*r1*(Lag2(5))*(1-R3)...
-s8*s11*(Lag2(9)+Lag2(13))*(1-R3);
%delay differential equation for the rate of change of the females between
%ages 8-11 R1 where Bf1(t-tau_2) is represented by Lag1(5), Cf1(t-tau_2)
%is represented by Lag1(9), Df1(t-tau_2) is represented by Lag1(13), Bf1(t-tau_3)
%is represented by Lag2(5), and Cf1(t-tau_3) is represented by Lag2(13).
dydt(6)=s8*(r1/2)*(Lag1(5)+Lag1(6))*Q3*Q4...
+s8*(Lag1(9)+Lag1(10)+Lag1(13)+Lag1(14))*Q3*Q4...
-a2*y(6)...
-s8*s11*(r1/2)*(Lag2(5)+Lag2(6))*R3*R4...
-s8*s11*(Lag2(9)+Lag2(10)+Lag2(13)+Lag2(14))*R3*R4;
%delay differential equation for the rate of change of females between ages
%8-11 R2 where Bf1(t-tau_2) is represented by Lag1(5), Bf2(t-tau_2) is
%represented by Lag1(6), Cf1(t-tau_2) is represented by Lag1(9),
%Cf2(t-tau_2) is represented by Lag1(10), Df1(t-tau_2) is represented by
%Lag1(13), Df2(t-tau_2) is represented by Lag1(14), Bf1(t-tau_3) is
%represented by Lag2(5), Bf2(t-tau_3) is represented by Lag2(6), Cf1(t-tau_3)
%is represented by Lag2(9), Cf2(t-tau_3) is represented by Lag2(10),
%Df1(t-tau_3) is represented by Lag2(13), Df2(t-tau_3) is represented by
%Lag2(14).
dydt(7)=s8*(r1/2)*(Lag1(5)+Lag1(6))*Q3*Q4...
+s8*(Lag1(9)+Lag1(10)+Lag1(13)+Lag1(14))*Q3*Q4...
-a2*y(7)...
-s8*s11*(r1/2)*(Lag2(5)+Lag2(6))*R3*R4...
-s8*s11*(Lag2(9)+Lag2(10)+Lag2(13)+Lag2(14))*R3*R4;
%delay differential equation for the rate of change of males between ages
%8-11 R2 where Bf1(t-tau_2) is represented by Lag1(5), Bf2(t-tau_2) is
%represented by Lag1(6), Cf1(t-tau_2) is represented by Lag1(9),
%Cf2(t-tau_2) is represented by Lag1(10), Df1(t-tau_2) is represented by
%Lag1(13), Df2(t-tau_2) is represented by Lag1(14), Bf1(t-tau_3) is
%represented by Lag2(5), Bf2(t-tau_3) is represented by Lag2(6), Cf1(t-tau_3)
%is represented by Lag2(9), Cf2(t-tau_3) is represented by Lag2(10),
%Df1(t-tau_3) is represented by Lag2(13), Df2(t-tau_3) is represented by
%Lag2(14).
dydt(8)=s8*Q6*r1*(Lag1(5)*Q3+Lag1(6))*(1-Q4)...
+s8*Q6*((Lag1(9)+Lag1(13))*Q3+Lag1(10)+Lag1(14))*(1-Q4)...
-a2*y(8)...
-s8*s11*R6*r1*(Lag2(5)*R3+Lag2(6))*(1-R4)...
-s8*s11*R6*((Lag2(9)+Lag2(13))*R3+Lag2(10)+Lag2(14))*(1-R4);
%delay differential equation for the rate of change of males between ages
%8-11 R3 where Bf1(t-tau_2) is represented by Lag1(5), Bf2(t-tau_2) is
%represented by Lag1(6), Cf1(t-tau_2) is represented by Lag1(9),
%Cf2(t-tau_2) is represented by Lag1(10), Df1(t-tau_2) is represented by
%Lag1(13), Df2(t-tau_2) is represented by Lag1(14), Bf1(t-tau_3) is
%represented by Lag2(5), Bf2(t-tau_3) is represented by Lag2(6), Cf1(t-tau_3)
%is represented by Lag2(9), Cf2(t-tau_3) is represented by Lag2(10),
%Df1(t-tau_3) is represented by Lag2(13), Df2(t-tau_3) is represented by
%Lag2(14).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 11-31 Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 11-31:y(9), R2 females between 11-31:y(10), R2 males between 11-31:y(11), R3 males between 11-31:y(12)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(9)=s8*s11*r1*(Lag2(5))*(1-R3)...
+s8*s11*(Lag2(9)+Lag2(13))*(1-R3)...
-a_new*y(9)...
-s8*s11*s31*r1*(Lag3(5))*(1-S3)...
-s8*s11*s31*(Lag3(9)+Lag3(13))*(1-S3);
%delay differential equation for the rate of change of females between ages
%11-31 R1 where Bf1(t-tau_3) is represented by Lag2(5), Cf1(t-tau_3) is
%represented by Lag2(9), Df1(t-tau_3) is represented by Lag2(13), Bf1(t-tau_5)
%is represented by Lag3(5), Cf1(t-tau_5) is represented by Lag3(9),
%Df1(t-tau_5) is represented by Lag3(13).
dydt(10)=s8*s11*(r1/2)*(Lag2(5)+Lag2(6))*R3*R4...
+s8*s11*(Lag2(9)+Lag2(10)+Lag2(13)+Lag2(14))*R3*R4...
-a_new*y(10)...
-s8*s11*s31*(r1/2)*(Lag3(5)+Lag3(6))*S3*S4...
-s8*s11*s31*(Lag3(9)+Lag3(10)+Lag3(13)+Lag3(14))*S3*S4;
%delay differential equation for the rate of change of females between ages
%11-31 R2 where Bf1(t-tau_3) is represented by Lag2(5), Bf2(t-tau_3) is
%represented by Lag2(6), Cf1(t-tau_3) is represented by Lag2(9),
%Cf2(t-tau_3) is represented by Lag2(10), Df1(t-tau_3) is represented by
%Lag2(13), Df2(t-tau_3) is represented by Lag2(14), Bf1(t-tau_5)
%is represented by Lag3(5), Bf2(t-tau_5) is represented by Lag3(6) Cf1(t-tau_5)
%is represented by Lag3(9), Cf2(t-tau_5) is represented by Lag3(10),
%Df1(t-tau_5) is represented by Lag3(13), and Df2(t-tau_5) is represented by Lag3(14).
dydt(11)=s8*s11*(r1/2)*(Lag2(5)+Lag2(6))*R3*R4...
+s8*s11*(Lag2(9)+Lag2(10)+Lag2(13)+Lag2(14))*R3*R4...
-a_new*y(11)...
-s8*s11*s31*(r1/2)*(Lag3(5)+Lag3(6))*S3*S4...
-s8*s11*s31*(Lag3(9)+Lag3(10)+Lag3(13)+Lag3(14))*S3*S4;
%delay differential equation for the rate of change of males between ages
%11-31 R2 where Bf1(t-tau_3) is represented by Lag2(5), Bf2(t-tau_3) is
%represented by Lag2(6), Cf1(t-tau_3) is represented by Lag2(9),
%Cf2(t-tau_3) is represented by Lag2(10), Df1(t-tau_3) is represented by
%Lag2(13), Df2(t-tau_3) is represented by Lag2(14), Bf1(t-tau_5)
%is represented by Lag3(5), Bf2(t-tau_5) is represented by Lag3(6) Cf1(t-tau_5)
%is represented by Lag3(9), Cf2(t-tau_5) is represented by Lag3(10),
%Df1(t-tau_5) is represented by Lag3(13), and Df2(t-tau_5) is represented
%by Lag3(14).
dydt(12)=s8*s11*R6*r1*(Lag2(5)*R3+Lag2(6))*(1-R4)...
+s8*s11*R6*((Lag2(9)+Lag2(13))*R3+Lag2(10)+Lag2(14))*(1-R4)...
-a_new*y(12)...
-s8*s11*s31*S6*r1*(Lag3(5)*S3+Lag3(6))*(1-S4)...
-s8*s11*s31*S6*((Lag3(9)+Lag3(13))*S3+Lag3(10)+Lag3(14))*(1-S4);
%delay differential equation for the rate of change of males between ages
%11-31 R3 where Bf1(t-tau_3) is represented by Lag2(5), Bf2(t-tau_3) is
%represented by Lag2(6), Cf1(t-tau_3) is represented by Lag2(9),
%Cf2(t-tau_3) is represented by Lag2(10), Df1(t-tau_3) is represented by
%Lag2(13), Df2(t-tau_3) is represented by Lag2(14), Bf1(t-tau_5)
%is represented by Lag3(5), Bf2(t-tau_5) is represented by Lag3(6) Cf1(t-tau_5)
%is represented by Lag3(9), Cf2(t-tau_5) is represented by Lag3(10),
%Df1(t-tau_5) is represented by Lag3(13), and Df2(t-tau_5) is represented
%by Lag3(14).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The 31+ Age Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 females between 11-31:y(9), R2 females between 11-31:y(10), R2 males between 11-31:y(11), R3 males between 11-31:y(12)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(13)=s8*s11*s31*r1*(Lag3(5))*(1-S3)...
+s8*s11*s31*(Lag3(9)+Lag3(13))*(1-S3)...
-a3*y(13);
%delay differential equation for the rate of change of females between ages
%31+ R1 where Bf1(t-tau_5) is represented by Lag3(5), Cf1(t-tau_5) is
%represented by Lag3(9), and Df1(t-tau_5) is represented by Lag3(13).
dydt(14)=s8*s11*s31*(r1/2)*(Lag3(5)+Lag3(6))*S3*S4...
+s8*s11*s31*(Lag3(9)+Lag3(10)+Lag3(13)+Lag3(14))*S3*S4...
-a3*y(14);
%delay differential equation for the rate of change of females between ages
%31+ R2 where Bf1(t-tau_5) is represented by Lag3(5), Bf2(t-tau_5) is
%represented by Lag3(6), Cf1(t-tau_5) is represented by Lag3(9),
%Cf2(t-tau_5) is represented by Lag3(10), Df1(t-tau_5) is represented by
%Lag3(13), and Df2(t-tau_5) is represented by Lag3(14).
dydt(15)=s8*s11*s31*(r1/2)*(Lag3(5)+Lag3(6))*S3*S4...
+s8*s11*s31*(Lag3(9)+Lag3(10)+Lag3(13)+Lag3(14))*S3*S4...
-a3*y(15);
%delay differential equation for the rate of change of males between ages
%31+ R2 where Bf1(t-tau_5) is represented by Lag3(5), Bf2(t-tau_5) is
%represented by Lag3(6), Cf1(t-tau_5) is represented by Lag3(9),
%Cf2(t-tau_5) is represented by Lag3(10), Df1(t-tau_5) is represented by
%Lag3(13), and Df2(t-tau_5) is represented by Lag3(14).
dydt(16)=s8*s11*s31*S6*r1*(Lag3(5)*S3+Lag3(6))*(1-S4)...
+s8*s11*s31*S6*((Lag3(9)+Lag3(13))*S3+Lag3(10)+Lag3(14))*(1-S4)...
-a3*y(16);
%delay differential equation for the rate of change of males between ages
%31+ R3 where Bf1(t-tau_5) is represented by Lag3(5), Bf2(t-tau_5) is
%represented by Lag3(6), Cf1(t-tau_5) is represented by Lag3(9),
%Cf2(t-tau_5) is represented by Lag3(10), Df1(t-tau_5) is represented by
%Lag3(13), and Df2(t-tau_5) is represented by Lag3(14).
dydt=dydt"; %the solver must return a column vector and this sets Matlab
%up to return a row vector, so we take the transpose to keep it
%all the same.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DDE23 requires histories for all variables even if they do not have a delay in
%the equation. Due to evidence the death die out after a certain amount of time we choose
%constants for our histories.
%This is a function file that maintains the system histories.
function s = Alli3Delay_hist_region(t) %function name and input
%for dde23, the output needs to be a column
%vector
%In this case we are assigning constant values to the variable histories.
hist_1 = 10; %this is the history of the females between ages 0-8 population from R1
hist_2 = 10; %this is the history of the females between ages 0-8 population from R2
hist_3 = 10; %this is the history of the males between ages 0-8 population from R2
hist_4 = 10; %this is the history of the males between ages 0-8 population from R3
hist_5 = 10; %this is the history of the females between ages 8-11 population from R1
hist_6 = 7; %this is the history of the females between ages 8-11 population from R2
hist_7 = 4; %this is the history of the males between ages 8-11 population from R2
hist_8 = 4; %this is the history of the males between ages 8-11 population from R3
hist_9=2.5; %this is the history of the females between ages 11-31 population from R1
hist_10 = 1.5; %this is the history of the females between ages 11-31 population from R2
hist_11=1; %this is the history of the males between ages 11-31 population from R2
hist_12=1; %this is the history of the males between ages 11-31 population from R3
hist_13=2.5; %this is the history of the females between ages 31+ population from R1
hist_14 = 1.5; %this is the history of the females between ages 31+ population from R2
hist_15=1; %this is the history of the males between ages 31+ population from R2
hist_16=1; %this is the history of the males between ages 31+ population from R3
s = [hist_1; hist_2; hist_3; hist_4; hist_5; hist_6; hist_7;hist_8;hist_9;
hist_10;hist_11;hist_12;hist_13;hist_14;hist_15;hist_16];
%This puts all the histories into a column vector and names it s. This is
%the format that dde23 would like it to be in.
end