dm 'log;clear;output;clear;odsresults;clear';

*******************************************************Create data******************************************;

data normal; 
seed =1976; 
call streaminit(seed); 
do obs = 1 to 501; 
e = rand ('normal',0,2); 
output; 
end; 
drop seed; 
run; 
 
data ma1; 
set normal; 
u=e+ 0.3*lag1(e);
run;



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
	drop obs e;
	If _N_ = 1 then delete;
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
