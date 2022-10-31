dm 'log;clear;output;clear;odsresults;clear';

*******************************************************Create data******************************************;

data normal; 
seed =1976; 
call streaminit(seed); 
do obs = 1 to 504; 
e0 = rand ('normal',0,2); 
output; 
end; 
drop seed; 
run; 

data error;
set normal;
e1=lag1(e0);
e2=lag2(e0);
run;

data normal1;
set error;
do i=1 to 504;
if _n_=1 then e=e0;
else if _n_=2 then e=e0+1.1*e1;
else e=e0+1.1*e1-0.28*e2;
end;
drop i e0 e1 e2 obs;
run;

proc transpose data=normal1 out=normal2;
run;

proc print data=normal2;run;

data ARMA; 
set normal2;
array e[*] col1-col504; 
array Y[504];
do i=1 to 504;   
j=i-1;
k=i-2;
if i=1 then Y[i]=e[i];
else if i=2 then Y[i]=0.6*Y[j]+e[i];
else Y[i]=0.6*Y[j]-0.25*Y[k]+e[i];
end;
drop col1-col504 Y1 Y2 i j k _NAME_;
run;

proc print data=ARMA;run;

proc transpose data=ARMA out=ARMA2;
run;

data ma1;
set ARMA2;
u=col1;
drop _NAME_ col1;
run;

proc print data=ma1;run;

*********************************************************Create lags ********************************************;
*Lag Macro*;
%macro lag(vars, lags);
	%let m = %sysfunc(countw(&vars));
	%do i=1 %to &m;
		%let var = %scan(&vars,&i);

		%do j=1 %to &lags;
			%do;
				lag_&var.&j = lag&j(&var);
			%end;
		%end;
	%end;
%mend lag;
*using lag macro to create lag columns;
data ma15;
	set ma1;
	%lag(u, 15);  /* variables to create lag_xxx to var xxx */
run;


*****************************************************Calculate Autocorrelations*********************************;
proc sql noprint; 
select mean(u) into :average from ma15; 
quit; 

data acf; 
set ma15;
array lags[*] lag_u1-lag_u15; 
array lagDif[15];
array lagx[15];
U0 = (u-&average); 
U02=U0**2;
do i=1 to 15;   
lagdif[i]=lags[i] - &average; 
lagx[i]=U0*lagdif[i];
end;
drop lagDif1-lagDif15 lag_u1-lag_u15 u U0 i;
run;

proc summary data=acf;
var lagx1--U02;
output out=totals sum=;
run;

data acf_res;
set totals;
array lags[*] lagx1-lagx15;
array acfs[15];
do i=1 to 15;   
acfs[i]=lags[i]/U02; 
end;
keep acfs1-acfs15;
run;

proc transpose data=acf_res out=want;
run;

data acfs; 
set want;
lags=_n_;
rename col1=ACF;
run;



******************************************Print Final Calculated ACF********************************;
proc print data=acfs;
var lags ACF;
run;


********************************************************CI calculation & Print CI with AFCs ********************************;
data acfs_cum;
  set acfs;
  p2=ACF**2;
  by lags notsorted;
   sumbylags+p2;
run;

data CI;
	set acfs_cum;
	 v=1/500*(1+p2);
	 sd=sqrt(v);
	 CIP=1.96*sd;
	 CIN=0-CIP;
run;

proc print data=CI;
var lags ACF CIP CIN;
run;


 *******************************************************CI Plot of AFCs*****************************************************;
proc sgplot data=CI;
band x=lags lower=CIN upper=CIP / transparency=0.5;
series x=lags y=ACF;
run;


proc sgplot data=CI;
  vbarparm category=lags response=ACF / barwidth=0.3 limitupper=CIP limitlower=CIN
           datalabel datalabelpos=top;
  yaxis grid;
run;


*******************************************************Using ready built code*****************************************************;
proc timeseries data=ma1 plots=corr outcorr=corr_pvals;
   var u; 
  corr lag n acf pacf acfstd acfprob pacfprob / nlag=15;
run;

proc print data=corr_pvals;
run;




*************************************************************Partial Autocorrelation****************************************;
*Noted slight difference with sas because sas uses  Yule-Walker approximation  method for PAC calculation;

 proc reg data=ma15 outest=est;
      m1: model u=lag_u1-lag_u15/noint;
run;
proc print data=est;
   run;

 Data pacf_res;
 set est;
 keep lag_u1-lag_u15;
 run;

proc transpose data=pacf_res out=want;
run;

data pacfs; 
set want;
lags=_n_;
rename col1=PACF;
run;



******************************************Print Final Calculated PACF********************************;
proc print data=pacfs;
var lags PACF;
run;


********************************************************CI calculation & Print CI with PAFCs ********************************;
data CI;
	set pacfs;
	 v=1/500;
	 sd=sqrt(v);
	 CIP=1.96*sd;
	 CIN=0-CIP;
run;

proc print data=CI;
var lags PACF CIP CIN;
run;


 *******************************************************CI Plot of AFCs*****************************************************;
proc sgplot data=CI;
band x=lags lower=CIN upper=CIP / transparency=0.5;
series x=lags y=PACF;
run;


proc sgplot data=CI;
  vbarparm category=lags response=PACF / barwidth=0.3 limitupper=CIP limitlower=CIN
           datalabel datalabelpos=top;
  yaxis grid;
run;
