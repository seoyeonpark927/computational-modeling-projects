(* ===== Double-Precision Computation ===== *)

ClearAll["Global`*"]
(*Define the function and experiment parameters*)
x0 = 0.5;
step[k_]:=10.^(-k)
f[x_]:=Exp[x]

(*Compute the approximations in double precision for each step size*)
doubleResults = Table[
Module[{h,xm, xc, xp, fm, f0, fp, a0, a1, a2, interpApprox, firstApprox, secondApprox},
h = step[k];

(* three nodes *)
xm = x0 - h;
xc = x0;
xp = x0 +h;

(* function samples *)
fm = f[xm];
f0 = f[xc];
fp = f[xp];

(* coefficients of the quadratic interpolant through the three nodes *)
a0 = f0;
a1 = 1/2 *(fp-fm);
a2 = 1/2*fm-f0+1/2*fp;

(* three approximations *)
interpApprox = a0+a1*0.5+a2*0.5^2;
firstApprox = a1/h;
secondApprox = 2*a2 / h^2;

{k, h, interpApprox, firstApprox, secondApprox}],
 {k,0,12}
];

doubleResults // TableForm

(*Extract step sizes and approximations from the double-precision results*)
hDouble = doubleResults[[All,2]];
interpDouble = doubleResults[[All,3]];
firstDouble = doubleResults[[All,4]];
secondDouble = doubleResults[[All,5]];

(*Compute the analytical reference values for error calculation*)
interpExact = Table[
Exp[0.5 + hDouble[[i]]/2], 
{i, 1, Length[hDouble]}
];
firstExact=Exp[0.5];
secondExact=Exp[0.5];

(*Compute and display the absolute errors in double precision*)
firstErrorDouble=Abs[firstDouble-firstExact];
secondErrorDouble=Abs[secondDouble-secondExact];
interpErrorDouble=Abs[interpDouble-interpExact];

{hDouble,interpErrorDouble,firstErrorDouble,secondErrorDouble}//TableForm

(*Organize the double-precision errors into a table*)
doubleErrors=Transpose[{hDouble,interpErrorDouble,firstErrorDouble,secondErrorDouble}];

doubleErrors//TableForm

(*Extract step sizes and interpolation errors*)
doubleErrors[[All,{1,2}]]

(* ===== Single-Precision Computation ===== *)

(*Convert values to simulated single precision (Real32)*)
to32[z_]:=First@Normal@NumericArray[{N[z]},"Real32"]

(*Define the function and experiment parameters in single precision*)
x0 = to32[0.5];
step[k_]:=to32[10.^(-k)]
f[x_]:=to32[Exp[x]]

(*Compute the approximations in single precision for each step size*)
singleResults = Table[
Module[{h,xm, xc, xp, fm, f0, fp, a0, a1, a2, interpApprox, firstApprox, secondApprox, term1, term2, sum1},
h = to32[step[k]];

(* three nodes *)
xm = to32[x0 - h];
xc = to32[x0];
xp = to32[x0 +h];

(* function samples *)
fm = to32[f[xm]];
f0 = to32[f[xc]];
fp = to32[f[xp]];

(* coefficients of the quadratic interpolant through the three nodes *)
a0 = to32[f0];
a1 = to32[to32[1/2] *(to32[fp-fm])];
a2 = to32[to32[to32[to32[1/2] *fm]-f0]+to32[to32[1/2] *fp]];

(* three approximations *)
term1=to32[a1*to32[0.5]];
	term2=to32[a2*to32[0.5^2]];
	sum1=to32[a0+term1];
	interpApprox=to32[sum1+term2];
firstApprox = to32[a1/h];
secondApprox = to32[to32[2*a2] / to32[h^2]];

{k, h, interpApprox, firstApprox, secondApprox}],
 {k,0,12}
];

singleResults//TableForm

(*Extract step sizes and approximations from the single-precision results*)
hSingle=singleResults[[All,2]];
interpSingle=singleResults[[All,3]];
firstSingle=singleResults[[All,4]];
secondSingle=singleResults[[All,5]];

(*Compute the analytical reference values for error calculation*)
interpExact = Table[
Exp[0.5 + hSingle[[i]]/2], 
{i, 1, Length[hSingle]}
];
firstExact=Exp[0.5];
secondExact=Exp[0.5];

(*Compute the absolute errors in single precision*)
firstErrorSingle=Abs[firstSingle-firstExact];
secondErrorSingle=Abs[secondSingle-secondExact];
interpErrorSingle=Abs[interpSingle-interpExact];

(*Organize the single-precision errors into a table*)
singleErrors=Transpose[
{hSingle,
interpErrorSingle,
firstErrorSingle,
secondErrorSingle}
];

singleErrors//TableForm

(* ===== Error Plots for Single and Double Precision ===== *)

(*Remove zero errors because they cannot be shown on a logarithmic scale*)
doubleInterpPlotData=Select[
doubleErrors[[All,{1,2}]],
#[[2]]>0& 
];
singleInterpPlotData=Select[singleErrors[[All,{1,2}]],#[[2]]>0&];
(*Plot interpolation errors for double and single precision*)
ListLogLogPlot[
{doubleInterpPlotData,singleInterpPlotData},
Joined->True,
PlotMarkers->Automatic,
PlotLegends->{"Double","Single"},
AxesLabel->{"h","Absolute Error"},
PlotLabel->"Interpolation Error"
]

(*Remove zero errors because they cannot be shown on a logarithmic scale*)
doubleFirstPlotData=Select[doubleErrors[[All,{1,3}]],#[[2]]>0&];
singleFirstPlotData=Select[singleErrors[[All,{1,3}]],#[[2]]>0&];
(*Plot first derivative errors for double and single precision*)
ListLogLogPlot[{doubleFirstPlotData,singleFirstPlotData},
Joined->True,
PlotMarkers->Automatic,
PlotLegends->{"Double","Single"},
AxesLabel->{"h","Absolute Error"},
PlotLabel->"First Derivative Error"
]

(*Remove zero errors because they cannot be shown on a logarithmic scale*)
doubleSecondPlotData=Select[doubleErrors[[All,{1,4}]],#[[2]]>0&];
singleSecondPlotData=Select[singleErrors[[All,{1,4}]],#[[2]]>0&];
(*Plot second derivative errors for double and single precision*)
ListLogLogPlot[{doubleSecondPlotData,singleSecondPlotData},
Joined->True,
PlotMarkers->Automatic,
PlotLegends->{"Double","Single"},
AxesLabel->{"h","Absolute Error"},
PlotLabel->"Second Derivative Error"
]

(* ===== Minimum Error Analysis ===== *)

MinimalBy[doubleErrors,#[[2]]&];
doubleInterpMin=MinimalBy[Select[doubleErrors,#[[2]]>0&],#[[2]]&];
doubleFirstMin=MinimalBy[Select[doubleErrors,#[[3]]>0&],#[[3]]&];
doubleSecondMin=MinimalBy[Select[doubleErrors,#[[4]]>0&],#[[4]]&];

doubleInterpMin
doubleFirstMin
doubleSecondMin

singleInterpMin=MinimalBy[Select[singleErrors,#[[2]]>0&],#[[2]]&];

singleFirstMin=MinimalBy[Select[singleErrors,#[[3]]>0&],#[[3]]&];

singleSecondMin=MinimalBy[Select[singleErrors,#[[4]]>0&],#[[4]]&];

singleInterpMin
singleFirstMin
singleSecondMin

(* ===== Floating-Point Node Collapse ===== *)

Table[
{hDouble[[i]],0.5+hDouble[[i]]==0.5,0.5-hDouble[[i]]==0.5},{i,1,Length[hDouble]}
]//TableForm

Table[
{hSingle[[i]],to32[to32[0.5]+hSingle[[i]]]==to32[0.5],to32[to32[0.5]-hSingle[[i]]]==to32[0.5]},{i,1,Length[hSingle]}
]//TableForm
