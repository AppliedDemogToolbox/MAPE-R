/*  There is no warranty for this code.
*/


/*  Enter the MAPER function into SAS
*/

%macro MAPER(data=, out=, var=, min=, max=, LAMBDA=);

/*  Usage:
      &DATA is the input dataset,
      &OUT is the output dataset to hold MAPE-R,
      &VAR is the variable for which to calculate MAPE-R,
      &MIN is the minimum value of LAMBDA to be used (-2 is recommended),
      &MAX is the maximum value of LAMBDA to be used (2 is recommended),
      &LAMBDA is the name of an optional dataset to hold the final value of LAMBDA.
*/

/*  Adapted from golden section search code from Coleman, Charles D., 
    "A Fast, High-Precision Implementation of the Univariate One-Parameter Box-Cox Transformation Using the Golden Section Search in SAS/IML®", 
    Proceedings of the 17th Northeast SAS Users Group (NESUG) Conference, 2004.  
    Available at http://www.nesug.org/proceedings/nesug04/an/an12.pdf.
*/

  proc iml;
    use &data;
	read all var {&var} into y;
	n=nrow(y);
	sumlog=sum(log(y));
	tol=1e-7;
    r=0.61803399;
    a=&min;
    b=&max;

	start llf(lambda) global(tol, y, n, sumlog);
      if abs(lambda)<tol then yl=log(y);
	  else yl=(y##lambda-1)/lambda;
	  avgyl=yl[:,];
	  f=-n*log(ssq(yl-avgyl)/n)/2+(lambda-1)*sumlog;
	  return(f);
	finish llf;
	
	fa=llf(a);
    fb=llf(b);
    c=a+r*(b-a);
    fc=llf(c);

	do while ((fc<fa) | (fc<fb));
	  if (fc<fa) then do;
	    a=a-2;
        fa=llf(a);
	  end;
	  else do;
	    b=b+2;
        fb=llf(b);
	  end;
	  c=a+r*(b-a);
      fc=llf(c);
    end;

	cdflag='c';
	do while (b-a>tol);
	  diff=r*(b-a);
	  if (cdflag='d') then do;
	    c=a+diff;
        fc=llf(c);
	  end;
      if (cdflag='c') then do;
        d=b-diff;
        fd=llf(d);
      end; 
      if (fc>fd) then do;
	    a=d;
        fa=fd;
        d=c;
        fd=fc;
        cdflag='d';
	  end;
      else do;
	    b=c;
        fb=fc;
        c=d;
        fc=fd;
        cdflag='c';
	  end;
    end;

    if (fa>fb) then lambda=a;
    else lambda=b;

    MAPER = (Y ## LAMBDA)[:,] ** (1/LAMBDA);

    %if &LAMBDA ne %then %do;
  	  create &lambda from lambda[colname='lambda'];
  	  append from lambda;
    %end;
    create &out from maper[colname='MAPE_R'];
    append from maper;
 
  quit;
%mend;


/*  Read data in 
    (example data is available at https://github.com/AppliedDemogToolbox/MAPE-R/raw/master/fictdata.xls)
*/

proc import DBMS=XLS
datafile="fictdata.xls" out=APES replace;
run;


/*  Apply the MAPER function to the data
*/

%MAPER(DATA=APES,OUT=TEST,VAR=APE,MIN=-2,MAX=2,LAMBDA=L);


/*  Print the output information
    (example output is available at https://github.com/AppliedDemogToolbox/MAPE-R/raw/master/mape-r.lst)
*/

proc print data=TEST;
proc print data=L;


