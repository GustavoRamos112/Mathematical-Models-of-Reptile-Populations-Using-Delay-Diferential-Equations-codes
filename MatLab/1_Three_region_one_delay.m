%RegionSolver1Delay_nd.m
%This is a script file that calls two function files, RegionDelayFun_nd,
%and RegionDelay_hist_nd. The purpose of this file is to numerically solve
%the 1 delay alligator model using Matlab"s built in dde solver, Runge
%Kutta Method, dde23.
%Comments updated for AMSSI 2007 on July 25, 2007


%%% 1. PARAMETER VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global a1 a2 c1 c2 b s  %declare variables as global to allow the dde23 solver
%to access their values, but we only have to change
%them in one place. The global command is to be
%used before you have defined the variables in the
%code.


b=.826;	%b is the birth rate and is set at a desired value
dJ=.2;	%dJ is the death rate of the juveniles and is set at a desired value 
dA=.0928;	%dA is the death rate of the adults and is set at a desired value
 
k1=.79;	%k1 is the carrying capacity of region 1 and is set at a desired value 
k2=.15;	%k2 is the carrying capacity of region 2 and is set at a desired value 
k3=.06;	%k3 is the carrying capacity of region 3 and is set at a desired value

a1=(dJ/b); %a1 is the death rate of juveniles divided by the birth rate 
a2=(dA/b); %a2 is the death rate of adults divided by the birth rate 
c1=(k2/k1); %c1 is the carrying capacity of region 2 divided by the carrying
%capacity of region 1
c2=(k3/k1); %c2 is the carrying capacity of region 3 divided by the carrying
%capacity of region 1
s=exp(-dJ*10);	%s is a probability of the survival to maturity



%%% 2. VARIABLE INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = 0;	%starting time value
tF = 200;	%the final time for which we will compute a solution
tspan = [t0 tF];%vector of initial and final time values for the dde23 solver.


tau=(10*b);	%set our time delay value to tau
lags=tau;	%"lags" is the array that holds the time delay values



%%% 3. NUMERICAL SOLUTION OF ODES %%%%%%%%%%%%%%%%%%%%%%%


Regionsol=dde23(@RegionDelayFun_nd,lags,@RegionDelay_hist_nd,tspan);
%this line calls the delay differential equation solver dde23.
%"RegionDelayFun_nd" is the name of the file storing our dde formulas, lags
%is the array of delay that we need to use, "RegionDelay_hist_nd" is the
%function file that maintains the system history and "tspan" feeds in our
%initial and final times for computation. The solutions to this system will
%be stored in the array"Regionsol".


%%% 4. PLOTTING OF SOLUTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%first, separate the output array into its various components


t = Regionsol.x;	%assigns the time outputs to a vector "t"
f1= Regionsol.y(1,:);	%assigns female juveniles in R1 outputs array
%to a vector named after the variable it
%represents, "f1".
f2= Regionsol.y(2,:);	%assigns female juveniles in R2 outputs array
%to a vector named after the variable it
 
%represents, "f2".

m2= Regionsol.y(3,:);	%assigns male juveniles in R2 outputs array
	%to a vector named after the variable it
	%represents, "m2".
m3= Regionsol.y(4,:);	%assigns male juveniles in R3 outputs array
	%to a vector named after the variable it
	%represents, "m3".
F1= Regionsol.y(5,:);	%assigns female adults in R1 outputs array
	%to a vector named after the variable it
	%represents, "F1"
F2= Regionsol.y(6,:);	%assigns female adults in R2 outputs array
	%to a vector named after the variable it
	%represents, "F2".
M2= Regionsol.y(7,:);	%assigns male adults in R2 outputs array
	%to a vector named after the variable it
	%represents, "M2".
M3= Regionsol.y(8,:);	%assigns male adults in R3 outputs array
	%to a vector named after the variable it
	%represents, "M3".



f=f1+f2+F1+F2;  %assigns the variable f to be equal to the total female
%population
m=m2+m3+M2+M3; %assigns the variable m to be equal to the total male
%population
r=(m./(m+f));	%assigns r to be equal to the total male population divided
%by the total male population plus the total female population
%will give you the ratio of males to total population


q=f1+F1;	%assigns q to be equal to the total population of R1
o=f2+F2;	%assigns o to be equal to the total population of females in R2 
k=m2+M2;	%assigns k to be equal to the total population of males in R2 
j=m3+M3;	%assigns j to be equal to the total population of R3



figure  %opens a new figure window

%subplot(2,1,1); %able to have 2 graphs in one window
hold on	%keeps all the following plots on the same graph
plot(t, f1, "b-", "LineWidth",3)	%plot the female juveniles in R1 a blue dashed line 
plot(t, f2, "g-", "LineWidth",3)	%plot the female juveniles in R2 a green dashed line 
plot(t, m2, "c-", "LineWidth",3)	%plot the male juveniles in R2 a light blue dashed line
plot(t, m3, "r-", "LineWidth",3)	%plot the male juveniles in R3 a red dashed line
%plot(t, F1, "b-", "LineWidth",3)	%plot the female adults in R1 a blue solid line
%plot(t, F2, "g-", "LineWidth",3)	%plot the female adults in R2 a green solid line
%plot(t, M2, "c-", "LineWidth",3)	%plot the male adults in R2 a light blue solid line
%plot(t, M3, "r-", "LineWidth",3)	%plot the male adults in R3 a red solid line 
xlabel("time")	% creates the label "time" on the x-axis 
ylabel("alligators") %creates the label "alligators" on the y axis
title("Three Regions One Delay") %creates a title for the graph 
legend("FR1","FR2","MR2","MR3")%"AFR1", "AFR2","AMR2","AMR3")	%creates a legend for the graph

figure  %opens a new figure window
hold on	%keeps all the following plots on the same graph
plot(t, q, "b-", "LineWidth",3)	%plot the total female population of R1 
plot(t, o, "g-", "LineWidth",3)	%plot the total female population of R2 
plot(t, k, "c-", "LineWidth",3)	%plot the total male population of R2 
plot(t, j, "r-", "LineWidth",3)	%plot the total male population of R3 
xlabel("time")	%creates the label "time" on the x-axis 
ylabel("alligators")	%creates the label "alligators" on the y-axis 
title("Two Delay Three Regions")	%creates a title for the graph 
legend("FR1", "FR2", "MR2", "MR3") %creates a legend for the graph

%subplot(2,1,2); %able to have 2 graphs in one window
figure
hold on	%this keeps all the following plots on the same graph 
plot(t, r, "m-", "LineWidth",3) %plots the ratio of males to total 
title("ratio of males to total population")	%creates the title
xlabel("time") %creates the label "time" on the x-axis
ylabel("alligators") %creates the label "alligators" on the y-axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a function code that stores all the delay differential equations
%formulas in another file to be called by the script file


function dydt=RegionDelayFun_nd(t,y,Z) %function name and variables, time
%variable and lags

global a1 a2 c1 c2 b s
Lag1=Z(:,1);	%a column vector whose rows represent the solution at the
%time delay for f1, f2, F1, F2, m2, m3, M2, and M3 respectively.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%The  Juvenile  Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 female juveniles:y(1), R2 female juveniles:y(2), R2 male juveniles:y(3), R3 male juveniles:y(4)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dydt(1) = (1/(1+y(5)))*y(5) - a1*y(1) - s*Lag1(5)*(1/(1+Lag1(5)));
%delay differential equation where F1(t-tau) is represented by Lag1(5).


dydt(2) = (1/2)*(c1/(c1+y(5)+y(6)))*(y(5)^2/(1+y(5)) + y(6))- a1*y(2) - (s/2)*((Lag1(5)^2/(1+Lag1(5)))+Lag1(6))*(c1/(c1+Lag1(5)+Lag1(6)));
%delay differential equation where F1(t-tau) is represented by Lag1(5) and
%F2(t-tau) is represented by Lag1(6).


dydt(3) = (1/2)*(c1/(c1+y(5)+y(6)))*(y(5)^2/(1+y(5)) + y(6))...
- a1*y(3) - (s/2)*((Lag1(5)^2/(1+Lag1(5)))+Lag1(6))*(c1/(c1+Lag1(5)+Lag1(6)));
%delay differential equation where F1(t-tau) is represented by Lag1(5) and
%F2(t-tau) is represented by Lag1(6).


dydt(4)  =  (c2/(c2+y(5)+y(6)))*(y(5)^2/(1+y(5))  +  y(6))*((y(5)+y(6))/(c1+y(5)+y(6)))...
- a1*y(4) - s*(c2/(c2+Lag1(5)+Lag1(6)))*((Lag1(5)^2)/(1+Lag1(5)) + Lag1(6))*((Lag1(5)+Lag1(6))/(c1+Lag1(5)+Lag1(6)));
%delay differential equation where F1(t-tau) is represented by Lag1(5) and
%F2(t-tau) is represented by Lag1(6)



%%%%%%%%%%%%%%%%%%%%%%The Adult Population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 female adults:y(5), R2 female adults:y(6), R2 male adults:y(7), R3 male adults:y(8)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dydt(5) = s*Lag1(5)*(1/(1+Lag1(5))) - a2*y(5); %delay differential equation, x"(t)%
%delay differential equation where F1(t-tau) is represented by Lag1(5)


dydt(6) = (s/2)*((Lag1(5)^2/(1+Lag1(5)))+Lag1(6))*(c1/(c1+Lag1(5)+Lag1(6))) - a2*y(6);
%delay differential equation where F1(t-tau) is represented by Lag1(5) and
%F2(t-tau) is represented by Lag1(6)


dydt(7) = (s/2)*((Lag1(5)^2/(1+Lag1(5)))+Lag1(6))*(c1/(c1+Lag1(5)+Lag1(6))) - a2*y(7);
%delay differential equation where F1(t-tau) is represented by Lag1(5) and
%F2(t-tau) is represented by Lag1(6)


dydt(8)  =  s*(c2/(c2+Lag1(5)+Lag1(6)))*((Lag1(5)^2)/(1+Lag1(5))  +  Lag1(6))*((Lag1(5)+Lag1(6))/(c1+Lag1(5)+Lag1(6)))...
- a2*y(8);
%delay differential equation where F1(t-tau) is represented by Lag1(5) and
%F2(t-tau) is represented by Lag1(6)
 
dydt=dydt'; %the solver must return a column vector and this sets Matlab
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


function s = RegionDelay_hist_nd(t) %function name and input
%for dde23, the output needs to be a column
%vector
%In this case we are assigning constant values to the variable histories.


f1_hist = 10;	%this is the history of the female juvenile population from R1 
f2_hist = 10;	%this is the history of the female juvenile population from R2 
m2_hist = 10;	%this is the history of the male juvenile population from R2 
m3_hist = 10;	%this is the history of the male juvenile population from R3 
F1_hist = 15;	%this is the history of the female adult population from R1 
F2_hist = 10;	%this is the history of the female adult population from R2 
M2_hist = 5;	%this is the history of the male adult population from R2 
M3_hist = 5;	%this is the history of the male adult population from R3



s = [f1_hist;f2_hist;m2_hist;m3_hist;F1_hist;F2_hist;M2_hist;M3_hist];

end
%This puts all the histories into a column vector and names it s. This is
%the format that dde23 would like it to be in.

