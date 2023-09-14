(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PatrickMLenggenhager`HyperBloch`"];


GetFullGraph;
GetSitePosition;
ImportCellGraphString;
ImportModelGraphString;


Begin["`Private`"];


(* ::Section:: *)
(*Initialization*)


Print["HyperBloch - Version 0.0.1\nMain author: Patrick M. Lenggenhager\n\nThis package loads the following dependencies:\n\t- L2Primitives by Srdjan Vukmirovic\n\t- NCAlgebra by J. William Helton and Mauricio de Oliveira"];


Get["L2Primitives.m"];


Quiet@Get["NCAlgebra`"];
SetOptions[inv, Distribute->True];
SetCommutative[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z];


(* ::Section:: *)
(*Definitions*)


(* ::Subsection:: *)
(*Helper Functions*)


GetFullGraph[graph_] := Graph[
	VertexList@graph,
	EdgeList@DirectedGraph[EdgeList@graph /. DirectedEdge -> UndirectedEdge],
	Sequence@@AbsoluteOptions[graph]
]


(* ::Subsection:: *)
(*Group Elements and Vertex Positions*)


InterpretGroupElementString[str_, rules_] := NCExpand@ToExpression@StringReplace[
	str,
	Join[{RegularExpression["\\^([\\-0-9]+)"]->"^($1)","*"->"**"}, (#[[1]] -> ToString[#[[2]]])&/@ rules]
];


(* Author: Tomas Bzdusek *)
Options[GetSitePosition] = {Orientation->"Default"};
GetSitePosition[tg_, fs_, expr_, OptionsPattern[]] := Module[
	{p, q, r, P, Q, R, ops, triangle, op, i, rules},
	(* triangle group signature *)
	{r, q, p} = tg;
	
	(* position of sites of fundamental triangle *)
	{P,Q,R} = Switch[OptionValue[Orientation],
		"Default",{
			PDPoint[{0,0}],
			PDPoint[\[Sqrt]((Cos[\[Pi] * (1/p+1/q)]+Cos[\[Pi]/r])/(Cos[\[Pi] * (1/p-1/q)]+Cos[\[Pi]/r])){1,0}],
			PDPoint[\[Sqrt]((Cos[\[Pi] * (1/p+1/r)]+Cos[\[Pi]/q])/(Cos[\[Pi] * (1/p-1/r)]+Cos[\[Pi]/q])){Cos[Pi/p],Sin[Pi/p]}]
		}
	];
	
	Block[{x, y, z, a, b, c},
		SetNonCommutative[x, y, z, a, b, c];
		rules = {"x" -> x, "y" -> y, "z" -> z, "a" -> a, "b" -> b, "c" -> c};
		((NCExpand[InterpretGroupElementString[expr, rules]^-1]) /. NonCommutativeMultiply -> Composition /. {
			1 -> (#&),
			x -> Function[{s}, L2Rotation[R, 2Pi/r][s]],
			y -> Function[{s}, L2Rotation[Q, 2Pi/q][s]],
			z -> Function[{s}, L2Rotation[P, 2Pi/p][s]],
			a -> Function[{s}, L2Reflection[LLine[{P, R}]][s]],
			b -> Function[{s}, L2Reflection[LLine[{Q, R}]][s]],
			c -> Function[{s}, L2Reflection[LLine[{P, Q}]][s]],
			Power[x, e_] :> Function[{s}, L2Rotation[R, 2Pi e/r][s]],
			Power[y, e_] :> Function[{s}, L2Rotation[Q, 2Pi e/q][s]],
			Power[z, e_] :> Function[{s}, L2Rotation[P, 2Pi e/p][s]],
			Power[a, e_] :> If[Mod[e, 2] == 1, Function[{s}, L2Reflection[LLine[{P,R}]][s]], (#&)],
			Power[b, e_] :> If[Mod[e, 2] == 1, Function[{s}, L2Reflection[LLine[{Q,R}]][s]], (#&)],
			Power[c, e_] :> If[Mod[e, 2] == 1, Function[{s}, L2Reflection[LLine[{P,Q}]][s]], (#&)]
		})[{R, Q, P}[[fs]]]
	]
]


(* ::Subsection:: *)
(*Import of Cell Graphs*)


ImportCellGraphString[str_, qname_]:=Module[{
		tg, specs, rels, center, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vlbls, vcoords,
		vertices, edges, etransls, faces, boundary
	},
	{tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vertices, edges, etransls, faces, boundary} =
		StringSplit[StringReplace[str,{"["->"{","]"->"}"}],"\n"];

	(* info *)
	tg = ToExpression@tg; (* triangle group signature *)
	{rels,center} = StringTrim[#, {"{", "}"}]&/@ StringSplit[specs, "},"];
	rels = StringSplit[rels, ", "]; (* cell relators *)
	center = ToExpression@center; (* cell center *)
	
	(* graph *)
	vertices = ToExpression@vertices;
	edges = DirectedEdge[vertices[[#1]], vertices[[#2]], #3]&@@@ToExpression[edges];
	faces = Table[
		Graph[
			Table[
				If[(edges[[face[[2,i]]]][[;;2]]/.DirectedEdge -> List) ==
						vertices[[face[[1, {i, Mod[i+1, Length[face[[1]]], 1]}]]]],
					edges[[face[[2,i]]]],
					edges[[face[[2,i]]]][[{2,1,3}]]
				],
				{i, 1, Length[face[[1]]]}
			]
		],
		{face, ToExpression[faces][[center]]}
	];
	
	(* algebra *)
	\[CapitalGamma]gens = "(" <> # <> ")" &/@(AssociationThread@@(StringTrim[StringSplit[StringTrim[#, {"{", "}"}], ","], " "]&/@StringSplit[\[CapitalGamma]gens, " -> "]));
	TD\[CapitalGamma] = StringTrim[StringSplit[StringTrim[TD\[CapitalGamma], {"{", "}"}], ","], " "];
	vlbls = (StringSplit[StringTrim[#, {"{ ", " }"}], ", "]&/@StringSplit[StringTrim[TGGw, {"{ ", " }"}], " }, { "])[[#[[1]], #[[2]]]]&/@vertices;
	
	(* translations *)
	etransls = StringTrim[StringSplit[StringTrim[etransls,{"{","}"}],","]," "];
	boundary = If[boundary=="{ }", {},
		{#1, #2, ToExpression@#3, ToExpression@#4, ToExpression@#5, #6}&@@@(
			StringTrim[StringSplit[#, ","], " "]&/@StringTrim[StringSplit[StringTrim[boundary, {"{", "}"}], "}, {"], {" { ", " } ", " "}]
		)
	];
	
	(* coordinates *)
	vcoords = Table[
		LToGraphics[GetSitePosition[tg, vertices[[i,1]], vlbls[[i]]], Model->PoincareDisk][[1]],
		{i, Length[vertices]}
	];
	
	<|
		"TriangleGroup" -> tg,
		"QuotientGroup" -> qname,
		"Graph"-> Graph[vertices, edges, VertexCoordinates -> vcoords],
		"VertexLabels" -> vlbls,
		"SchwarzTriangleLabels" -> TD\[CapitalGamma],
		"EdgeTranslations" -> etransls,
		"BoundaryEdges" -> boundary,
		"TranslationGenerators" -> \[CapitalGamma]gens,
		"Faces" -> faces
	|>
]


(* ::Subsection:: *)
(*ImportModelGraph*)


ImportModelGraphString[str_, qname_]:=Module[{
		tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, model, vertices, edges, etransls, faces,
		rels, center, vlbls, vcoords
	},
	{tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, model, vertices, edges, etransls, faces} =
		StringSplit[StringReplace[str, {"["->"{", "]"->"}"}], "\n"];
	
	(* info *)
	tg = ToExpression@tg; (* triangle group signature *)
	{rels, center} = StringTrim[#, {"{", "}"}]&/@ StringSplit[specs, "},"];
	rels = StringSplit[rels, ", "]; (* cell relators *)
	center = ToExpression@center; (* cell center *)
	
	(* algebra *)
	\[CapitalGamma]gens = "(" <> # <> ")" &/@(AssociationThread@@(StringTrim[StringSplit[StringTrim[#, {"{", "}"}], ","], " "]&/@StringSplit[\[CapitalGamma]gens, " -> "]));
	TD\[CapitalGamma] = StringTrim[StringSplit[StringTrim[TD\[CapitalGamma], {"{", "}"}], ","], " "];
	
	(* graph *)
	vertices = ToExpression@vertices;
	edges = DirectedEdge[vertices[[#1]], vertices[[#2]], #3]&@@@ToExpression[edges];
	faces = Table[
		Graph[
			Table[
				If[(edges[[face[[2,i]]]][[;;2]]/.DirectedEdge -> List) ==
						vertices[[face[[1, {i, Mod[i+1, Length[face[[1]]], 1]}]]]],
					edges[[face[[2,i]]]],
					edges[[face[[2,i]]]][[{2,1,3}]]
				],
				{i, 1, Length[face[[1]]]}
			]
		],
		{face, ToExpression[faces]}
	];
	
	(* vertex labels and coordinates *)
	vlbls = (StringSplit[StringTrim[#, {"{ ", " }"}], ", "]&/@StringSplit[StringTrim[TGGw, {"{ ", " }"}], " }, { "])[[#[[1]], #[[2]]]]&/@vertices;	
	vcoords = Table[
		LToGraphics[GetSitePosition[tg, vertices[[i,1]], vlbls[[i]]], Model->PoincareDisk][[1]],
		{i, Length[vertices]}
	];
	
	(* edge translations *)
	etransls = StringTrim[StringSplit[StringTrim[etransls,{"{", "}"}], ","], " "];
	
	<|
		"TriangleGroup" -> tg,
		"QuotientGroup" -> qname,
		"Genus" -> ToExpression@StringReplace[qname,RegularExpression["T(\\d+)\\.(\\d+)"]->"$1"],
		"Graph" -> Graph[vertices, edges, VertexCoordinates -> vcoords],
		"UndirectedGraph"->UndirectedGraph@Graph[vertices, edges, VertexCoordinates -> vcoords],
		"FullGraph" -> GetFullGraph@Graph[vertices, edges, VertexCoordinates -> vcoords],
		"VertexLabels" -> vlbls,
		"SchwarzTriangleLabels" -> TD\[CapitalGamma],
		"EdgeTranslations" -> etransls,
		"TranslationGenerators" -> \[CapitalGamma]gens,
		"Faces" -> faces
	|>
]


(* ::Section:: *)
(*Package Footer*)


End[];
EndPackage[];
