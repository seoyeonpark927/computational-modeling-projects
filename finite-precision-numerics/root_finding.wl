(* ===== Bisection method for F(x) = x^3 - x - 1 on [1, 2] ===== *)

(*Define the function and initialize the bisection method*)
ClearAll["Global`*"]

F[x_]:=x^3-x-1
a=1.;
b=2.;
count=0;
bisHistory={};

{a,b,count,F[a],F[b]}

(*Perform bisection iterations until the desired tolerance is reached*)
While[b-a>2*10^-10,
m = (a+b)/2;
AppendTo[bisHistory,m];
If[F[a]*F[m]<0,
b = m,
a=m
];
count = count + 1;
]

{a, b, count, F[a], F[b]}

(*Compute and display the final bisection results*)
xApprox=(a+b)/2;
finalBracket={a,b};
residual=Abs[F[xApprox]];
xExact=1.32471795724474602596;
absoluteError=Abs[xApprox-xExact];

{
{"Final bracket",finalBracket},
{"Approximation",xApprox},
{"Residual",residual},
{"Absolute error",absoluteError},
{"Bracket width",b-a},
{"Number of halvings",count}
}//TableForm

(*Display the final bracket and approximation with full machine precision*)

InputForm[a]

InputForm[b]

InputForm[xApprox]

(* ===== Fixed-point iteration x -> (1 + x)^(1/3)  (contraction) ===== *)

(*Define the fixed-point function and initialize the iteration*)
g[x_]:=(1+x)^(1/3);

xOld=1.;
xNew=g[xOld];
countFP=1;
fpHistory = {xOld, xNew};

(*Perform fixed-point iterations until the desired tolerance is reached*)
While[
Abs[xNew-xOld]>10^-10&&countFP<200,
xOld = xNew;
xNew = g[xOld];
AppendTo[fpHistory,xNew];
countFP=countFP+1;
]

(*Display the final fixed-point approximation and iteration count*)

InputForm[xNew]

Abs[xNew-xOld]

countFP

xApproxFP=xNew;

(*Compute and display the final fixed-point results*)

residualFP=Abs[F[xApproxFP]];
xExact=1.32471795724474602596;
absoluteErrorFP=Abs[xApproxFP-xExact];

{
{"Approximation",xApproxFP},
{"Residual",residualFP},
{"Absolute error",absoluteErrorFP},
{"Successive difference",Abs[xNew-xOld]},
{"Iterations",countFP}
}//TableForm

InputForm[xApproxFP]

(*Compute the absolute errors at each iteration*)
fpErrors=Abs[fpHistory-xExact];
fpErrors

(*Compare the convergence of bisection and fixed-point iteration*)

fpPlotData=Transpose[{Range[0,Length[fpErrors]-1],fpErrors}];
fpPlotData

bisErrors=Abs[bisHistory-xExact];
bisPlotData=Transpose[{Range[1,Length[bisErrors]],bisErrors}];

ListLogPlot[
{bisPlotData,fpPlotData},
Joined->True,
PlotMarkers->Automatic,
PlotLegends->{"Bisection","Fixed-Point Iteration"},
AxesLabel->{"Iteration","Absolute Error"},
PlotLabel->"Convergence Comparison",
GridLines->Automatic,
ImageSize->Large
]

(* ===== Fixed-point iteration y -> y^3 - 1  (not a contraction: diverges) ===== *)

(*Define the alternative fixed-point function and initialize the iteration*)
ClearAll[yOld,yNew,yHistory,countC,g2];

g2[y_]:=y^3-1;

yOld=1.;
countC=0;
yHistory={yOld};

(*Perform the alternative iteration and record its behavior*)
While[
countC<20&&Abs[yOld]<=10^6,
yNew = g2[yOld];
AppendTo[yHistory, yNew];

yOld = yNew;
countC=countC+1;
]

yHistory
countC
InputForm[yOld]
