Print["hello"]

(*(*inputDir="/Users/iryna/Documents/isotonic/SharedBennoIryna/CodingHTA2/20160727Data_2009Classif/"*)*)
inputDir=$CommandLine[[4]]
Print[inputDir];
(* ::Input:: *)
(*(*second argument to be passed:*)*)


(* ::Input:: *)
(*(*fToBinarize="20160727_varSortedThresh0.5_GeneOnly_2653_67528_15SDvs27O_2009classif_42p.txt"*)*)
fToBinarize= $CommandLine[[5]]
Print[fToBinarize];

(* ::Section:: *)
(*Functions*)


readStandardExcelMatrix[data_,opts___]:={First/@Rest[data],Rest[First[data]],Rest/@Rest[data]};


(*Clear[writeStandardExcelMatrix];
writeStandardExcelMatrix[{columns_,rows_,data_},filename_]:=Export[filename,Prepend[Transpose[Prepend[Transpose[data],rows]],Prepend[columns,0]],"CSV"] *)


(*PAR[opts___]:=(Print[opts];opts);*)


diffuse[data_,dx_]:=Flatten/@Transpose/@({#-Abs[#]*dx,#+Abs[#]*dx}&/@data)


WritePermutationData[data_,dx_,filename_]:=Module[{n,str,perm,theWords,words,perWord,perms},
perms=Ordering[Ordering[#]]&/@data;
(* If a row of data is (x_1, x_2, ... , x_n) *)
(* "perms" is (i_1, i_2, ..., i_n), where x_{i_1} <= x_{i_2} <= ... <= x_{i_n} *)
str=OpenWrite[filename,BinaryFormat->True];

(*Write header*)
BinaryWrite[str,n=Length[perms[[1]]],"UnsignedInteger64"];
BinaryWrite[str,Length[perms],"UnsignedInteger64"];

(*Write unmodified permutations*)
perWord=Floor[64/Log[2,n]];
words=Ceiling[n/perWord];
(BinaryWrite[str,(FromDigits[#,n]&/@Reverse/@(ArrayReshape[#-1,{words,perWord}])),"UnsignedInteger64"])&/@perms;
Print["Wrote ",Length[perms]," permutations."];

(*Write diffused permutations*)
perms=Ordering[Ordering[#]]&/@diffuse[data,dx];
n= 2n;
perWord=Floor[64/Log[2,n]];
words=Ceiling[n/perWord];
(BinaryWrite[str,(FromDigits[#, n]&/@Reverse/@(ArrayReshape[#-1,{words,perWord}])),"UnsignedInteger64"])&/@perms;
Print["Wrote ",Length[perms]," permutations."];
Close[filename];
];


(*ReadPermutationData[filename_]:=Module[{n,perWord,words,result1, theFormat, str, result2,numPerms},
str=OpenRead[filename,BinaryFormat\[Rule]True];
n=BinaryRead[str,theFormat="UnsignedInteger64"];
numPerms=BinaryRead[str,theFormat="UnsignedInteger64"];
perWord=Floor[64/Log[2,n]];
words=Ceiling[n/perWord];
Print[{n,numPerms,perWord,words}];
result1=ArrayReshape[BinaryReadList[str,theFormat,words*numPerms],{numPerms,words}];
result1=Take[(Join@@(Reverse[IntegerDigits[#,n,perWord]]&/@#)+1),n]&/@result1;
Print["Read ",Length[result1]," permutations."];

n=2n;
perWord=Floor[64/Log[2,n]];
words=Ceiling[n/perWord];
Print[{n,numPerms,perWord,words}];
result2=ArrayReshape[BinaryReadList[str,theFormat,words*numPerms],{numPerms,words}];
Print[result2[[1,1]]];
result2=Take[(Join@@(Reverse[IntegerDigits[#,n,perWord]]&/@#)+1),n]&/@result2;
Print["Read ",Length[result2]," permutations."];

Close[str];
{result1,result2}
];
*)


(* ::Input:: *)
(*(*WritePhenotypes[phenotypes_,filename_]:=Module[{n,str},*)
(*str=OpenWrite[filename];*)
(*Write[str,Length[phenotypes]];*)
(*Write[str,#]&/@phenotypes;*)
(*Print["Wrote ",Length[phenotypes]," phenotypes."];*)
(*Close[filename];*)
(*];*)*)



(*Action!*)

lFiles=FileNames[fToBinarize,inputDir]
Print[lFiles];
({exRows,exColumns,exData}=readStandardExcelMatrix[Import[#,"TSV"]];WritePermutationData[exData,0.03,#<>".bin"])&/@lFiles
