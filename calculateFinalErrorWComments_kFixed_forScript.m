(* ::Package:: *)

(* ::Section:: *)
(*Set folder and file names, other parameters *)


(* ::Text:: *)
(*Determine user so we can set right base path*)


(*iAmIryna=($UserName\[NotEqual]"benno")*)


(* ::Text:: *)
(*The variable describing which analysis project we're working on*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*projectSubfolderName="20160825AllViruses_last_longList/";*)*)
(**)


(* ::Text:: *)
(*Set base folder where inputs and cross-validation results can be found*)


(* ::Input:: *)
(*(*mainBaseFolder=If[iAmIryna,*)
(*"/Users/iryna/Documents/isotonic/SharedBennoIryna/DREAM/",*)
(*"/Users/benno/BTSync/iryna-benno/DREAM/"];*)*)
(**)


(* ::Input:: *)
(**)


projectBaseFolder=$CommandLine[[4]]
Print["hey1"]
Print[projectBaseFolder]


(* ::Input:: *)
(**)


bestK=ToExpression[$CommandLine[[5]]]
Print["hey2"]
Print[bestK]



(* ::Input:: *)
(**)


maxVotes=bestK;


(* ::Text:: *)
(*Names for subfolder containing concatenated and sorted files from cluster action/internal CV/CV/full data*)


(* ::Input:: *)
(**)


inputSubfolder="data/";
crossValidationSubfolder="isotonicResults/";
concatenatedAndSortedSubfolder="_concatenatedAndSorted/";
internalCVSubfolderPrefix="cvInternalSets"; 
cvSubfolderPrefix="cvTestSet";
realSubfolderPrefix="real";


(* ::Text:: *)
(*Set inputs folder name*)


(* ::Input:: *)
(*(*projectBaseFolder=mainBaseFolder<>projectSubfolderName*)*)


(* ::Input:: *)
(**)


inputBaseFolder=projectBaseFolder<>inputSubfolder;
Print[inputBaseFolder]


(* ::Subsubsection:: *)
(*Transcriptome data*)


transcriptsFile=FileNames["*EXPRESSION*",inputBaseFolder][[1]];
Print["hey3"]


(* ::Subsubsection:: *)
(*Annotations (needs to be fixed)*)


(*annotationsFile=FileNames["*ANNOTATIONS*",inputBaseFolder][[1]]*)


(* ::Subsubsection:: *)
(*Phenotypes*)


phenotypeFile=FileNames["*PHENOTYPE*",inputBaseFolder][[1]];


(* ::Subsubsection:: *)
(*Cross-validation results base folder*)


(* ::Input:: *)
(**)


crossValidationResultsBaseFolder=projectBaseFolder<>crossValidationSubfolder;


(* ::Subsubsection:: *)
(*Internal cross-validation results (leave 1-3-out)*)


(* ::Input:: *)
(*internalCVResultFolder=crossValidationResultsBaseFolder<>internalCVSubfolderPrefix<>concatenatedAndSortedSubfolder*)


(* ::Input:: *)
(*internalCVResultPattern="*Fromisoton*LOInt*.txt"*)


(* ::Subsubsection:: *)
(*Global CV cross-validation results (leave-1-out)*)


(* ::Input:: *)
(**)


cvResultFolder=crossValidationResultsBaseFolder<>cvSubfolderPrefix<>concatenatedAndSortedSubfolder;


(* ::Input:: *)
(**)


cvResultPatternSuffix="*LOTest_once.txt";


(* ::Subsubsection:: *)
(*Global result features*)
(**)


(* ::Input:: *)
(**)


globalResultsFolder=crossValidationResultsBaseFolder<>realSubfolderPrefix<>concatenatedAndSortedSubfolder;


(* ::Input:: *)
(**)


globalResultsPattern="Fromisoton_all_once.txt";


(* ::Subsubsection:: *)
(*Simple-leave-out-features*)


(* ::Input:: *)
(*simpleLeaveOutFolder=crossValidationResultsBaseFolder<>realSubfolderPrefix<>"/"<>internalCVSubfolderPrefix<>concatenatedAndSortedSubfolder*)


(* ::Section:: *)
(*Functions*)


(* ::Input:: *)
(**)


numberOfLinesInFile[f_String]:=Length[Import[f,"Lines"]];


readStandardExcelMatrix[data_,opts___]:={First/@Rest[data],Take[First[data],-Length[Rest[data[[2]]]]],Rest/@Rest[data]};


pointsX[{i_,j_,signI_,signJ_}]:=Transpose[{exData[[i]]*signI,exData[[j]]*signJ,phenotypes}]


pointsXCvCount[exDataCvCount_,phenotypesCvCount_,{i_,j_,signI_,signJ_}]:=Transpose[{exDataCvCount[[i]]*signI,exDataCvCount[[j]]*signJ,phenotypesCvCount}]


(* ::Subsection::Closed:: *)
(*Isotonic regression v6*)


isotonicRegression[allpointsInput_List,fun_:Less,fitOnly_List:{}]:=Block[
{i,c,cA,iMax,x,y,n,arrayInd,g,sums,columnToHs,cg,iParent,regSort,l,rightMin,reg,delta,minAttained,arrayToH,h,rightMinNode,yl,ylPlus1,taskStack={},theseIndices,lowReg,highReg,lowSet,highSet,regValue,allRegSort,thePerm,allReg,allpoints,rotDelta,hs,allPointPerm,lowInd,highInd,lastInLowSet,lastInHighSet,oldHScore,update,sib,change,ignoreInds,ignoreFlags,notIgnoreInds,changes,altScore,points,inputOrder,result},
allpoints=Sort[allpointsInput];
inputOrder=Ordering[Ordering[allpointsInput]];
ignoreInds=inputOrder[[fitOnly]];
pe["ignoreInds"];
ignoreFlags=Table[False,{Length[allpoints]}];
(ignoreFlags[[#]]=True)&/@ignoreInds;
pe["ignoreFlags"];
notIgnoreInds=Complement[Range@Length@allpoints,ignoreInds];
allReg=Last/@allpoints;
allRegSort=Sort[Union[allReg[[notIgnoreInds]]]];
regValue=Table[0,{Length[allpoints]}];(* output *)
allPointPerm=Range[Length[allpoints]];
pe["allPointPerm"];
altScore=0;
AppendTo[taskStack,{{1,Length[allpoints]},{1,Length[allRegSort]}}];
While[taskStack=!={},
{{{lowInd,highInd},{lowReg,highReg}},taskStack}={First[taskStack],Rest[taskStack]};
theseIndices=allPointPerm[[lowInd;;highInd]];
pe["allPointPerm"];
pe["Length[theseIndices]"];pe["theseIndices"];
If[theseIndices=={},Continue[]];
If[(lowReg==highReg),(regValue[[#]]=allRegSort[[lowReg]])&/@theseIndices; 
(altScore+=regValue[[#]]-allRegSort[[lowReg]])&/@theseIndices;
Continue[]];
points=Drop[#,-1]&/@allpoints[[theseIndices]];
reg=allReg[[theseIndices]];
regSort=Take[allRegSort,{lowReg,highReg}]; 
n=c=Length[points];
columnToHs=Ordering[Ordering[Last/@points]];
c=Length[columnToHs];
cA=2c+1;
arrayInd=RotateRight[Range[c+1]+c,rotDelta=BitShiftLeft[BitClear[c,BitLength[c]-1]+1]];
arrayToH=Join[Range[c],RotateLeft[Range[c+1],rotDelta]];
sums=Table[0,{cA}];
l=BitShiftRight[lowReg+highReg];
yl=allRegSort[[l]];ylPlus1=allRegSort[[l+1]];
pe["{yl,ylPlus1,allRegSort}"];
minAttained=Table[0,{n}];
pe["columnToHs"];
pe["columnToHs.Range[Length@columnToHs]"];
For[g=1,g<=n,g++,

(*First pass: calculate Z(g-1,h)*) 
cg=columnToHs[[g]];
pe["cg"];
rightMin=oldHScore=0; 
i=rightMinNode=arrayInd[[cg]]; 
update=Abs[reg[[g]]-{yl,ylPlus1}]; (*right, left*)
pe["arrayInd[[cg]]"];
pe["theseIndices[[g]]"];
pe["ignoreFlags[[theseIndices[[g]]]]"];
If[TrueQ@ignoreFlags[[theseIndices[[g]]]],update={0,0}];
pe["update"];
While[i>1,
oldHScore+=sums[[i]];
rightMin+=sums[[i]];
If[(BitAnd[i,1]==0) && (fun[sums[[sib=i+1]],rightMin]),
rightMinNode=sib;rightMin=sums[[sib]]];
i=BitShiftRight[i];
];

(*Second pass: Update tree for new score of g*)
change=update[[2]]+rightMin-oldHScore;
pe["change"];
changes={};
i=arrayInd[[cg]];
While[i>1,
sib=BitXor[i,1];
sums[[i]]+=change;
sums[[sib]]+=update[[BitAnd[i,1]+1]];
change=Min[sums[[i]],sums[[sib]]];
AppendTo[changes,{update[[BitAnd[i,1]+1]],change,i}];
sums[[i]]-=change;
sums[[sib]]-=change;
i=BitShiftRight[i];
];
pe["changes"];
pe["sums.Range[Length[sums]]"];

(*Third pass: Find min. k's by going down from rightMinNode *) 
i=rightMinNode;
While[i<=c,
i=BitShiftLeft[i];
If[fun[sums[[i+1]],sums[[i]]],i++];
];
minAttained[[g]]=i;
PrintMaybe["Sums: ",Plus@@Drop[sums,1]];
pe["altScore"];
];

(* Extract total minimum*) i=1;While[i<=c,i=BitShiftLeft[i];If[fun[sums[[i+1]],sums[[i]]],i++]];
PrintMaybe["Total minimum node: ", i];
lowSet={};highSet={};
g=n+1;h=arrayToH[[i]];
lastInLowSet=lowInd-1;lastInHighSet=highInd+1;
While[(--g>0),
If[h<=columnToHs[[g]],
--lastInHighSet;
PrependTo[highSet,g],

++lastInLowSet;
PrependTo[lowSet,g]
]; 
If[h==columnToHs[[g]],h=arrayToH[[minAttained[[g]] ]] ];
];
pe["theseIndices[[lowSet]]"];
pe["theseIndices[[highSet]]"];
allPointPerm=Join[allPointPerm[[;;lowInd-1]],theseIndices[[lowSet]],theseIndices[[highSet]],allPointPerm[[highInd+1;;]]];
(* Put two subproblems on task stack *)
PrependTo[taskStack,{{lowInd,lastInLowSet},{lowReg,l}}];
PrependTo[taskStack,{{lastInHighSet,highInd},{l+1,highReg}}];
];

(* Output regression solution *)
result=replaceLastCoordinate[allpoints,regValue,ignoreInds];
{First[result][[inputOrder]],Last[result]}
]


Clear[replaceLastCoordinate];
replaceLastCoordinate[inputPoints_List,fittedValues_List,ignoreInds_List:{}]:=Module[{totalDifference},
(* Compute total difference between input and fitted *)
totalDifference=Abs[(Last/@inputPoints)-fittedValues];
{Transpose[Append[Drop[Transpose[inputPoints],-1] ,fittedValues]],{Plus@@(totalDifference[[ignoreInds]]),Plus@@totalDifference}}
];


isotonicRegression[l1_List,l2_List,values_List]:=Block[{},
isotonicRegression[Transpose[{l1,l2,values}]]
];


(* ::Subsection::Closed:: *)
(*Randomization function for isotonicRegression input data*)


isotonicRegressionRandomizePositions[allpointsInput_List,xvalPos_List:{}]:=Block[{newData=allpointsInput},
(newData[[#,-1]]=RandomInteger[{0,1}])&/@xvalPos;
newData];


(* ::Subsection::Closed:: *)
(*Settings and functions*)


classesToValues={"Asympto"->0,"DF"->1,"DHF1"->8,"DSS"->9}


valueRange={0,1};meanValue=Mean[valueRange];


colorMap=ColorData[{"TemperatureMap",valueRange}];


readStandardExcelMatrix[data_,opts___]:={First/@Rest[data],Rest[First[data]],Rest/@Rest[data]};


{turquoise=RGBColor[0.,0.53,0.94],myYellow=RGBColor[0.91,0.9,0],myOrange=RGBColor[1.,0.6,0.23]}


subsample[l_List, r_]:=RandomSample[l,Floor[Length[l]*r]];


Clear[plotSituation];
plotSituation[points_List,opts___]:=ListPlot[Map[Drop[#,-1]&,SortBy[GatherBy[points,Last],#[[1,-1]]&],{2}],PlotLegends->theLegend,opts]
(*theClasses,AspectRatio\[Rule]1,PlotStyle\[Rule]{turquoise,myYellow,myOrange,Red}*)


signedMult[l_List,dx_]:=If[First[l]>0,l*dx,l/dx];


diffusePoints[points_List,dx_]:=Module[{firsts,seconds,vals},
{firsts,seconds,vals}=Transpose[points];
Join[Transpose[{signedMult[firsts,dx],signedMult[seconds,1/dx],vals}],Transpose[{signedMult[firsts,1/dx],signedMult[seconds,dx],vals}]]
]


theLegend=PointLegend[colorMap/@{0,0,0,3},theClasses,LegendMarkerSize->20]
(*theLegend=PointLegend[colorMap/@{0,1,2,3},theClasses,LegendMarkerSize\[Rule]20]*)



(*theLegend=PointLegend[colorMap/@{0,3},{"Non-DSS","DSS"},LegendMarkerSize\[Rule]20]*)


plotSituationContinuous[points_List,colorMap_,opts___]:=ListPlot[List/@(Most/@points),PlotStyle->colorMap/@(Last/@points),PlotLegends->theLegend,opts]


borderVisualization[inputPoints_List,colorMap_]:=Module[{fittedPoints,score,minmax,p1,p2},
{fittedPoints,score}=isotonicRegression[inputPoints];
minmax=Plus@@isotonicMinMax[fittedPoints];
p1=Graphics[backgroundPlot[fittedPoints,colorMap]];
p2=plotSituationContinuous[inputPoints,colorMap];
Print["Score: ",score];
Show[p1,p2]
];


(* ::Section:: *)
(*Read phenotypes, set numPatients*)


phenotypes=ToExpression[First[#]]&/@Rest[Import[phenotypeFile,"TSV"]];


(* ::Input:: *)
(**)


numPatients=Length[phenotypes] 


(* ::Section::Closed:: *)
(*Read annotations*)


(*
Dimensions[dataR=Import[annotationsFile,"TSV"]]


geneColumn=Last@First@Position[dataR,"PRTN3"]


annotRange={geneColumn-2,geneColumn}


{annotRows,annotColumns,orgAnnotData}=readStandardExcelMatrix[dataR];


Table[annotIndex[annotRows[[i]]]=i,{i,Length[annotRows]}];
*)



(* ::Input:: *)
(*orgAnnotData*)


(* ::Section:: *)
(*Read expression data*)


dataR=Import[transcriptsFile,"TSV"];


Dimensions/@({exRows,exColumns,exData}=readStandardExcelMatrix[dataR]);


(*Dimensions[annotData=orgAnnotData[[Table[annotIndex[exRows[[i]]],{i,Length[exRows]}]]]]*)


suffix="";


(exIndex[exColumns[[#]]]=#)&/@Range[Length[exColumns]];


Clear[class];
(class[First[#]]=Last[#])&/@Transpose[{exColumns,phenotypes}];


(* ::Section:: *)
(*More functions that use annotations/expression data*)


(* ::Subsection:: *)
(*Exploratory functions*)


Clear[isotonicMinMax];
isotonicMinMax[points_,keptPoints_]:=Module[{mins,maxs,orderedPoints,xOrd,yOrder=Ordering[points[[All,2]]],
n=Length[points],p,rumin,ip,newRow,rumax,min,max,newValue,isFitted,test},
isFitted=Table[True,{n}];
(isFitted[[#]]=False)&/@(InversePermutation[yOrder])[[keptPoints]];
{min,max}={Min[#],Max[#]}&[Last/@points[[keptPoints]]];
rumin={Table[min,{n}]};
For[i=1,i<=n,i++,
newValue=If[isFitted[[i]],min,Last[points[[yOrder[[i]]]]]];
newRow=PadRight [Table[min,{yOrder[[i]]-1}],n,newValue];
AppendTo[rumin,Max/@Transpose[{newRow,Last[rumin]}]];
];
rumin=Rest[rumin];
rumax={Table[max,{n}]};
For[i=n,i>0,i--,
newValue=If[isFitted[[i]],max,Last[points[[yOrder[[i]]]]]];
newRow=PadLeft[Table[max,{n-yOrder[[i]]}],n,newValue];
PrependTo[rumax,Min/@Transpose[{newRow,First[rumax]}]];
];
rumax=Most[rumax];

{Map[Flatten,Transpose[{Transpose@Outer[{##}&,First/@points,Sort[points[[All,2]]]],rumin},{3,2,1}],{2}],
Map[Flatten,Transpose[{Transpose@Outer[{##}&,First/@points,Sort[points[[All,2]]]],rumax},{3,2,1}],{2}]}
];


Clear[backgroundPlotOld];
backgroundPlotOld[points_,colorMap_,fittedPoints_]:=Module[{cornerValues,keptPoints,n=Length[points],i,j,shortPoints},
keptPoints=Complement[Range[n],fittedPoints];
cornerValues=Plus@@(isotonicMinMax[points,keptPoints]/2);
shortPoints=cornerValues[[All,All,1;;2]];
Flatten[#,2]&[Join[Table[{Polygon[{shortPoints[[i,j]],shortPoints[[i+1,j]],shortPoints[[i+1,j+1]],shortPoints[[i,j+1]]},VertexColors->((colorMap[Last[#]])&/@{cornerValues[[i,j]],cornerValues[[i+1,j]],cornerValues[[i+1,j+1]],cornerValues[[i,j+1]]})]},{i,n-1},{j,n-1}],
{Black,Circle[#,{0.2,0.2}]}&/@points[[keptPoints,1;;2]],{Green,Thickness[0.005],Circle[#,{0.05,0.05}]}&/@points[[fittedPoints,1;;2]]]]
];


Clear[backgroundPlot];
backgroundPlot[points_,colorMap_,fittedPoints_]:=Module[{cornerValues,keptPoints,n=Length[points],i,j,shortPoints},
keptPoints=Complement[Range[n],fittedPoints];
cornerValues=Plus@@(isotonicMinMax[points,keptPoints]/2);
shortPoints=cornerValues[[All,All,1;;2]];
Print[cornerValues[[1,1]]];
Flatten[#,2]&[Join[Table[{Polygon[{shortPoints[[i,j]],shortPoints[[i+1,j]],shortPoints[[i+1,j+1]],shortPoints[[i,j+1]]},VertexColors->((colorMap[Last[#]])&/@(If[Length@Union[Last/@#]==1,#,Table[{meanValue},{4}]]&@{cornerValues[[i,j]],cornerValues[[i+1,j]],cornerValues[[i+1,j+1]],cornerValues[[i,j+1]]}))]},{i,n-1},{j,n-1}],
{Black,Circle[#,{0.2,0.2}]}&/@points[[keptPoints,1;;2]],{Green,Thickness[0.005],Circle[#,{0.05,0.05}]}&/@points[[fittedPoints,1;;2]]]]
];


Clear[borderVisualizationNew];
borderVisualizationNew[inputPointsR_List,colorMap_,crossVal_:{0,1},opts___]:=Module[{fittedPoints,fittedPoints1,score,minmax,p1,p2,cvN,cvK,toFitPoints,n=Length[inputPointsR],keptPoints,inputPoints=Sort[inputPointsR],xscore},
{cvK,cvN}=crossVal;
Table[
toFitPoints=RandomSample[Range[n],cvK];
(*TEMP toFitPoints=Range[cvK];*)
fittedPoints=First@isotonicRegression[inputPoints,Less,toFitPoints];
(*fittedPoints+=First@isotonicRegression[inputPointsR,LessEqual,toFitPoints];
fittedPoints/=2;*)
{xscore,score}=Last@replaceLastCoordinate[inputPoints,Last/@fittedPoints,toFitPoints];
p1=Graphics[backgroundPlot[fittedPoints,colorMap,toFitPoints]];
p2=plotSituationContinuous[inputPointsR,colorMap];
{{xscore,score},{p1,p2}}
,{cvN}]
];


correctDengueValues[points_List]:={#[[1]],#[[2]],#[[3]]-If[#[[3]]>2,6,0]}&/@points;


Clear[isoReport];
isoReport[inputPointsR_,{xValK_Integer},statsQ_:False]:=isoReport[inputPointsR,Subsets[Range[Length[inputPointsR]],{xValK}],statsQ];isoReport[inputPointsR_,{xValK_Integer,xValN_Integer},statsQ_:False]:=isoReport[inputPointsR,Table[RandomSample[Range[Length@inputPointsR],xValK],{xValN}],statsQ];isoReport[inputPointsR_,concreteXVals_List,statsQ_:False]:=Module[{inputPoints=Sort@inputPointsR,xval,stats,n=Length@inputPointsR,xvalMean,globalMean,reg1,reg2,fittedPoints,p1,p2,xvalTotal,globalTotal},
stats=Last[isotonicRegression[inputPoints,Less,#]]&/@concreteXVals;
{xvalMean,globalMean}=Mean/@Transpose[N@stats];
{xvalTotal,globalTotal}=Total/@Transpose[N@stats];
If[!statsQ,
Print["Total cross validation error over ",Length[concreteXVals], " cross validations is ",xvalTotal," (Mean= ",xvalMean,")."];
Print["Total global fitting error is ",globalTotal,"(mean:",globalMean,")."];
];
If[statsQ,Return[{xvalMean,globalMean}]];
reg1=isotonicRegression[inputPoints,Less];
reg2=isotonicRegression[inputPoints,LessEqual];
(*fittedPoints=(First[reg1]+First[reg2])/2;*)
fittedPoints=First[reg1];
p1=Graphics@backgroundPlot[correctDengueValues[fittedPoints],colorMap,{}];
p2=plotSituationContinuous[correctDengueValues@inputPointsR,colorMap];
{p1,p2}
];


pointsReport[r_List,s_List,xvals_List]:=Module[{},
isoReport[pointsRX[r,s,phenotypes],xvals]
];

(*
pairReport[{a_,b_,signA_,signB_},xvals_List,statsQ_:False]:=Module[{isoRep},
If[!statsQ,
If [a<=Length[annotData],Print["Feature 1: ",TableForm[{Append[Take[annotData[[a]],annotRange],exRows[[a]]]}]]];
If [b<=Length[annotData],Print["Feature 2: ",TableForm[{Append[Take[annotData[[b]],annotRange],exRows[[b]]]}]]];
];
isoRep=isoReport[pointsX[{a,b,signA,signB}],xvals,statsQ];
If[!statsQ,Return[isoRep]];
Join[{a,b,signA,signB},isoRep]
];


pairReportDiffused[{a_,b_,signA_,signB_},xvals_List,dx_]:=Module[{},
If [a<=Length[annotData],Print["Feature 1: ",TableForm[{Append[Take[annotData[[a]],annotRange],exRows[[a]]]}]]];
If [b<=Length[annotData],Print["Feature 2: ",TableForm[{Append[Take[annotData[[b]],annotRange],exRows[[b]]]}]]];
isoReport[diffusePoints[pointsX[{a,b,signA,signB}],dx],xvals]
];


Last/@correctDengueValues[pointsIncludingHDLP[40939]]


ToExpression[StringSplit["2544 11884 1 -1"]]


plotOpts={Axes->True,PlotLegends->theClasses}


*)


(* ::Section:: *)
(*Compute average error for all cv patients and globally (function definitions)*)


(* ::Text:: *)
(*Define the size of the moving average interval (which is halfMovingAverageSize*2 + 1)*)


(* ::Input:: *)
(**)


halfMovingAverageSize=2;


(* ::Text:: *)
(*Set function to use for isotonic regression (Less or LessEqual)*)


(* ::Input:: *)
(**)


isotonicCompareFunction=Less;


(* ::Text:: *)
(*Set default value for maxVotes (will be reduced automatically as needed)*)


(* ::Input:: *)
(*(*maxVotes=100000;*)*)


(* ::Text:: *)
(*Determine global cv error for different numbers of voters, and performances when using best k determined systematically from internal leave-outs*)


(* ::Input:: *)
(**)


 


(* ::Input:: *)
(**)


calculateFinalError:=
Module[{allFileNames,filenameProfiles,matchingLengthFileNames,lengthMatchingPositions,goodPositions,cvPatient,patternToMatch,theFileNames,phenotypesCvCount,exDataCvCount,exColumnsCvCount,numCases,numCrossValidations,xvalSets,voterSets,individualVotes,ensembleErrors,numVoters,movingAverages,cvVoters,loocvAndEnsembleErrors,cvSetSize,j},

loocvAndEnsembleErrors=Transpose[Table[ 
patternToMatch="*isoton_"<>ToString[cvPatient]<>cvResultPatternSuffix;
theFileNames=Select[FileNames[patternToMatch,cvResultFolder],(numberOfLinesInFile[#]>0)&];
Assert[Length[theFileNames]==1];
cvVoters=(Drop[ToExpression[StringSplit[#]],2]&/@Take[First/@Import[theFileNames[[1]],"TSV"],bestK]);
maxVotes=Min[Min[numberOfLinesInFile/@theFileNames],maxVotes];
individualVotes=Table[(Last/@(First@isotonicRegression[pointsX[cvVoters[[i]]+{1,1,0,0}],isotonicCompareFunction,{cvPatient+1}]))[[cvPatient+1]],{i,bestK}];
{Abs[Round[Mean[individualVotes]]-phenotypes[[cvPatient+1]]],Mean[individualVotes]},
{cvPatient,0,numPatients-1}]];

Assert[Length[theFileNames=FileNames[globalResultsPattern,globalResultsFolder]]==1];
cvVoters=Take[Drop[ToExpression[StringSplit[#]],2]&/@(First/@Import[theFileNames[[1]],"TSV"]),maxVotes];
individualVotes=ParallelTable[(Last/@(First@isotonicRegression[pointsX[cvVoters[[i]]+{1,1,0,0}],isotonicCompareFunction,{j}]))[[j]],{j,numPatients},{i,maxVotes}];
{ensembleErrors=N@ParallelTable[Mean@Table[Mean[(Abs[(Round[Table[individualVotes[[j,i]],{i,numVoters}]]-phenotypes[[j]])])],{j,numPatients}],{numVoters,1,maxVotes}],loocvAndEnsembleErrors}
];

 
 
 
 
 
 
 



(* ::Input:: *)
(*calculateSimpleLeaveOutError:=*)
(*Module[{patternToMatch,theFileNames,allFileNames,cvSetSize,matchingLengthFileNames,goodPositions,filenameProfiles,phenotypesCvCount,numCrossValidations,xvalSets,voterSets,individualVotes,ensembleErrors,numVoters,bestK,movingAverages},*)
(*(* Redundant statement [SORRY :-<]; theFileNames is assigned again below *)*)
(*theFileNames=FileNames[internalCVResultPattern,simpleLeaveOutFolder];*)
(*(* Get file names that match pattern from directory *)*)
(*allFileNames=allFileNames=Select[FileNames[internalCVResultPattern,internalCVResultFolder],(numberOfLinesInFile[#]>0)&];*)
(*(* Get file profiles (as described in "calculateFinalError") *)*)
(*filenameProfiles=(Select[ToExpression[Rest@StringSplit[#,"_"]],NumberQ])&/@allFileNames;*)
(*(* Select size we're interested in (as described in "calculateFinalError") *)*)
(*cvSetSize=First[First[Tally[Length/@(filenameProfiles)]]];*)
(*(* Determine "lengthMatchingPositions" = indices of files with matching profiles size (as described in "calculateFinalError") *)*)
(*(* Same inefficiency as in "calculateFinalError" ["matchingLengthFileNames" not really needed] *)*)
(*matchingLengthFileNames=allFileNames[[lengthMatchingPositions=First/@Position[Length/@filenameProfiles,cvSetSize]]];*)
(*(* Determine only those profiles that match value of cvPatient  *)*)
(*goodPositions=lengthMatchingPositions[[First/@Position[First/@filenameProfiles[[lengthMatchingPositions]],cvPatient]]];*)
(*(* Determine actual filenames of interest now  *)*)
(*theFileNames=allFileNames[[goodPositions]];*)
(*(* It goes on just as "calculateFinalError" from here *)*)
(*xvalSets=Rest/@filenameProfiles[[goodPositions]]+1;*)
(*numCrossValidations=Length[theFileNames];*)
(*maxVotes=Min[Min[numberOfLinesInFile/@theFileNames],maxVotes];*)
(*voterSets=Table[(Drop[ToExpression[StringSplit[#]],2]&/@Take[First/@Import[theFileNames[[j]],"TSV"],maxVotes]),{j,numCrossValidations}];*)
(*individualVotes=ParallelTable[(Last/@(First@isotonicRegression[pointsX[voterSets[[j,i]]+{1,1,0,0}],isotonicCompareFunction,xvalSets[[j]]]))[[xvalSets[[j]]]],{j,numCrossValidations},{i,maxVotes}];*)
(*ensembleErrors=ParallelTable[Mean@Table[Mean[(Abs[(Round[Mean[Table[individualVotes[[j,i]],{i,numVoters}]]]-phenotypes[[xvalSets[[j]]]])])],{j,numCrossValidations}],{numVoters,maxVotes}];*)
(*bestK=First[Ordering[ movingAverages=MovingAverage[ensembleErrors,Table[1,{halfMovingAverageSize*2+1}]],1]]+halfMovingAverageSize;*)
(*Print["Best k=",bestK, " (error=",N[movingAverages[[bestK-halfMovingAverageSize]],3],") derived from ",numCrossValidations, " cross-validations."];*)
(*movingAverages*)
(*];*)


(* ::Input:: *)
(**)


(* ::Subchapter:: *)
(*DREAM*)


(* ::Section:: *)
(*calculateFinalError: Compute average error for all cv patients*)


(* ::Text:: *)
(*Define the size of the moving average interval (which is halfMovingAverageSize*2 + 1)*)


(* ::Input:: *)
(**)


halfMovingAverageSize=2;


(* ::Text:: *)
(*Set the total number of patients left out in internal cross-validation (this determines which filenames we're looking at from which we parse the left-out patients)*)


(* ::Input:: *)
(**)


isotonicCompareFunction=Less;
Print["Starting To Calculate Serious Stuff"]


(* ::Input:: *)
(**)


result=calculateFinalError;


(* ::Input:: *)
(**)


Print["Ended Serious Stuff"]
phenotypePredictions=result[[2]][[2]];
(**)


(* ::Input:: *)
(**)


result[[1]];


(* ::Input:: *)
(**)


Export[projectBaseFolder<>"phenotypeProbabilities.tsv", N[phenotypePredictions], "TSV"];






(* ::Input:: *)
(*(**)
(**)
(*(*I'm generating the predictions for new patients*)*)
(*patientsToPredict={1,2};*)
(*bestK=7;(* get individual votes when we leave out only patient j *)*)
(*(*theFileNames=FileNames[globalResultsFolder, globalResultsPattern]*)*)
(*theFileNames=List[globalResultsFolder<>globalResultsPattern]*)
(*(* get the first maxVotes voters *)*)
(*cvVoters=Take[Drop[ToExpression[StringSplit[#]],2]&/@(First/@Import[theFileNames[[1]],"TSV"]),maxVotes]*)
(*individualVotes=ParallelTable[(Last/@(First@isotonicRegression[pointsX[cvVoters[[i]]+{1,1,0,0}],isotonicCompareFunction,{j}]))[[j]],{j,numPatients},{i,bestK}]*)
(*PredictedPhenotypesLeaderBoard=Mean[Transpose[individualVotes]]*)*)
