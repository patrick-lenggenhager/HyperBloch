(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PatrickMLenggenhager`HyperBloch`"];


GetFullGraph;
GetSitePosition;

ImportCellGraphString;
ImportModelGraphString;

GetTriangleTessellation;
ShowTriangles;

GetSchwarzTriangle;
GetVertex;
GetEdge;
GetCellGraphVertex;
GetCellGraphEdge;
GetTranslatedCellGraphEdge;
GetCellGraphFace;

GetCellBoundary;

ShowCellGraph;
ShowCellSchwarzTriangles;
ShowCellGraphFlattened;
ShowCellBoundary;

VisualizeCellGraph;


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


(* ::Subsection::Closed:: *)
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


(* ::Subsection::Closed:: *)
(*Import of Cell Graphs*)


ImportCellGraphString[str_, qname_]:=Module[{
		tg, specs, rels, center, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vlbls, vcoords,
		vertices, edges, etransls, facesstr, faces,faceedges, boundary
	},
	{tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vertices, edges, etransls, facesstr, boundary} =
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
				If[e[[2]] == 1,
					edges[[e[[1]]]],
					edges[[e[[1]]]][[{2, 1, 3}]]
				],
				{e, face}
			]
		],
		{c, 1, 3}, {face, ToExpression[facesstr][[c]]}
	];
	faceedges =Map[edges[[#[[1]]]]&,ToExpression[facesstr][[center]], {2}];
	
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
		"CellCenter" -> center,
		"Graph"-> Graph[vertices, edges, VertexCoordinates -> vcoords],
		"VertexLabels" -> vlbls,
		"SchwarzTriangleLabels" -> TD\[CapitalGamma],
		"EdgeTranslations" -> etransls,
		"BoundaryEdges" -> boundary,
		"TranslationGenerators" -> \[CapitalGamma]gens,
		"Faces" -> faces[[center]],
		"FaceEdges" -> faceedges,
		"AllFaces" -> faces
	|>
]


(* ::Subsection::Closed:: *)
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
				If[e[[2]] == 1,
					edges[[e[[1]]]],
					edges[[e[[1]]]][[{2, 1, 3}]]
				],
				{e, face}
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


(* ::Subsection:: *)
(*Graphical Visualization*)


(* ::Subsubsection::Closed:: *)
(*Triangle Tessellations*)


(* Tomas' code *)
ClearAll[GetTriangleTessellation];
Options[GetTriangleTessellation]={ColorBoundary->RGBColor[{0,0,0,0}],ColorFill->LightGray,LineThickness->0.00125};
(* Specify (p,q,r) of the triangle group \[CapitalDelta](p,q,r). *)
(* The p-fold symmetric point is placed at the center, and *)
(* the q-fold symmetric point is placed on the horizontal axis. *)
GetTriangleTessellation[{r_,q_,p_},gens_,opts:OptionsPattern[]]:=GetTriangleTessellation[{r,q,p},gens,opts]=Module[{
(* Specify how many times to iterate the triangle generation. *)
(* Note that the number of triangles grows exponentially with this parameter! *)
Generations=gens,
(* Specify coloring of the triangle tessellation *)
ColorBoundary=OptionValue[ColorBoundary],
ColorFill=OptionValue[ColorFill],
LineThickness=OptionValue[LineThickness],

P,Q,R,tr,TriangleLength,TriangleList,cc,temp0,aa,temp1,bb,
temp2,TempLen,boundaries,ToFill,triangles,absolut
},

(* Some analytic results for the initial triangle. *)
P=PDPoint[{0,0}];
Q=PDPoint[Sqrt[(Cos[\[Pi] (1/p+1/q)]+Cos[\[Pi]/r])/(Cos[\[Pi] (1/p-1/q)]+Cos[\[Pi]/r])]*{1,0}];
R=PDPoint[Sqrt[(Cos[\[Pi] (1/p+1/r)]+Cos[\[Pi]/q])/(Cos[\[Pi] (1/p-1/r)]+Cos[\[Pi]/q])]*{Cos[Pi/p],Sin[Pi/p]}];



(* Specify first triangle. *)
tr=LLine[{P,Q,R,P}];
TriangleLength={0};
TriangleList={tr};



(* First iteration of the algorithm. *)
cc=1;

temp0=TriangleList[[1]];
For[aa=1,aa<p,aa++,
temp1=Chop@L2Rotation[P,aa*(2*Pi/p)][temp0];
TriangleList=AppendTo[TriangleList,temp1];
];

TriangleLength=AppendTo[TriangleLength,Length[TriangleList]];

(* All subsequent iterations of the algorithm. *)
For[cc=2,cc<Generations+1,cc++,

For[bb=TriangleLength[[cc-1]]+1,bb<TriangleLength[[cc]]+1,bb++,
temp1=TriangleList[[bb]];
For[aa=1,aa<p,aa++,
temp2=Chop@L2Rotation[temp1[[1,1]],aa*(2*Pi/p)][temp1];
TriangleList=AppendTo[TriangleList,temp2]
];
For[aa=1,aa<q,aa++,
temp2=Chop@L2Rotation[temp1[[1,2]],aa*(2*Pi/q)][temp1];
TriangleList=AppendTo[TriangleList,temp2]
];
For[aa=1,aa<r,aa++,
temp2=Chop@L2Rotation[temp1[[1,3]],aa*(2*Pi/r)][temp1];
TriangleList=AppendTo[TriangleList,temp2]
];
];

TempLen=Length[TriangleList];
For[aa=TempLen,aa>TriangleLength[[1]],aa--,
For[bb=1,bb<aa,bb++,
If[TriangleList[[aa]]==TriangleList[[bb]],{
TriangleList=Delete[TriangleList,aa];
Break[];
}
];
];
];

TriangleLength=AppendTo[TriangleLength,Length[TriangleList]];
];



(* All triangles have been generated. The remainder of *)
(* the code has to do with simple data manipulation and plotting. *)

boundaries = Graphics[{ColorBoundary,Thickness[LineThickness],LToGraphics[{Flatten@TriangleList}, Model->PoincareDisk]}];

ToFill=Table[LPolygon[{
TriangleList[[aa,1,1]],
TriangleList[[aa,1,2]],
TriangleList[[aa,1,3]]
}],{aa,1,TriangleLength[[-1]]}];

triangles= Graphics[{ColorFill,LToGraphics[{Flatten@ToFill}, Model->PoincareDisk]}];

Show[{triangles,boundaries}, PlotRange->{{-1,1}, {-1,1}}, AspectRatio->Automatic]
]


ClearAll[ShowTriangles]
Options[ShowTriangles] = {
	RasterizeGraphics -> True,
	NumberOfGenerations -> 2
};
ShowTriangles[tg_, opts:OptionsPattern[{ShowTriangles, Graphics, Rasterize, GetTriangleTessellation}]] := ShowTriangles[tg, opts] = Graphics[
	{
		If[OptionValue[RasterizeGraphics],
			Inset[Rasterize[
					Show[
						GetTriangleTessellation[tg, OptionValue[NumberOfGenerations],
							Sequence@@FilterRules[{opts}, Options[GetTriangleTessellation]]
						],
						Sequence@@FilterRules[{opts}, {ImageSize, PlotRange}],
						PlotRange -> 1.1{{-1, 1}, {-1, 1}},
						ImageSize -> 500
					],
					Sequence@@FilterRules[{opts}, Options[Rasterize]],
					ImageSize -> 500
				], Scaled[{0.5,0.5}], Scaled[{0.5,0.5}], Scaled[1]
			],
			GetTriangleTessellation[tg, OptionValue[NumberOfGenerations],
				Sequence@@FilterRules[{opts}, Options[GetTriangleTessellation]]
			][[1]]
		],
		Circle[]
	},
	Sequence@@FilterRules[{opts}, Options[Graphics]],
	PlotRange -> 1.1{{-1, 1}, {-1, 1}},
	ImageSize -> 500
]


(* ::Subsubsection::Closed:: *)
(*Cell Graph Elements*)


GetSchwarzTriangle[tg_, g_] := LToGraphics[
	LPolygon[GetSitePosition[tg, #, g]&/@{1, 2, 3}],
	Model -> PoincareDisk
]
GetVertex[tg_, v_] := LToGraphics[
	GetSitePosition[tg, v[[1]], v[[2]]],
	Model -> PoincareDisk
]
GetEdge[tg_, e_] := LToGraphics[
	LLine[GetSitePosition[tg, #[[1]], #[[2]]]&/@e],
	Model -> PoincareDisk
]


ResolveVertex[cgraph_, vertex_] := {
	vertex[[1]],
	cgraph["VertexLabels"][[Position[VertexList@cgraph["Graph"], vertex][[1, 1]]]]
}


ResolveTranslation[transl_, gens_] := StringReplace[transl,
	RegularExpression["g(\\d+)"] :> gens["g$1"]
]


ResolveEdge[cgraph_, edge_] := Module[{\[Gamma], v1, v2},
	If[MemberQ[EdgeList@cgraph["Graph"], edge],
		(* edge with default orientation *)
		\[Gamma] = ResolveTranslation[
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge][[1, 1]]]],
			cgraph["TranslationGenerators"]
		];
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];,
		(* edge with inverted orientation *)
		\[Gamma] = "(" <> ResolveTranslation[
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge[[{2, 1, 3}]]][[1, 1]]]],
			cgraph["TranslationGenerators"]
		] <> ")^(-1)";
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];
	];
	
	{v1, {v2[[1]], v2[[2]]<>"*"<>\[Gamma]}}
]


GetCellGraphVertex[cgraph_, vertex_] := GetVertex[
	cgraph["TriangleGroup"],
	ResolveVertex[cgraph, vertex]
]
GetCellGraphEdge[cgraph_, edge_] := GetEdge[
	cgraph["TriangleGroup"],
	ResolveEdge[cgraph, edge]
]


ResolveTranslatedEdge[cgraph_, edge_, \[Gamma]0_] := Module[{\[Gamma], v1, v2},
	If[MemberQ[EdgeList@cgraph["Graph"], edge],
		(* edge with default orientation *)
		\[Gamma] = ResolveTranslation[
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge][[1, 1]]]],
			cgraph["TranslationGenerators"]
		];
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];,
		(* edge with inverted orientation *)
		\[Gamma] = "(" <> ResolveTranslation[
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge[[{2, 1, 3}]]][[1, 1]]]],
			cgraph["TranslationGenerators"]
		] <> ")^(-1)";
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];
	];
	
	{{{v1[[1]], v1[[2]]<>"*"<>\[Gamma]0}, {v2[[1]], v2[[2]]<>"*"<>\[Gamma]<>"*"<>\[Gamma]0}}, \[Gamma]<>"*"<>\[Gamma]0}
]
GetTranslatedCellGraphEdge[cgraph_, edge_, \[Gamma]_] := GetEdge[
	cgraph["TriangleGroup"],
	ResolveTranslatedEdge[cgraph, edge, ResolveTranslation[\[Gamma], cgraph["TranslationGenerators"]]][[1]]
]


ResolveFace[cgraph_, face_] := Module[{\[Gamma] = "1", re},
	Table[
		re = ResolveTranslatedEdge[cgraph, e, \[Gamma]];
		\[Gamma] = re[[2]];
		re[[1]],
		{e, EdgeList[face]}
	]
]
GetCellGraphFace[cgraph_, face_] :=
LToGraphics[
	LPolygon[GetSitePosition[cgraph["TriangleGroup"], #[[1,1]], #[[1,2]]]&/@ResolveFace[cgraph, face]],
	Model -> PoincareDisk
]
(*GetEdge[cgraph["TriangleGroup"], #]&/@ResolveFace[cgraph, face]*)


(* ::Subsubsection::Closed:: *)
(*Cell Boundary*)


GetCellBoundary[cgraph_] := {
	#[[4]],
	LToGraphics[LLine[{
		GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[1, 1]], #[[1]]],
		GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[2, 1]], #[[2]]]
	}], Model -> PoincareDisk],
	LToGraphics[GetSitePosition[cgraph["TriangleGroup"], 3,
		ResolveTranslation[#[[6]], cgraph["TranslationGenerators"]]
	], Model -> PoincareDisk],
	#[[6]]
}&/@cgraph["BoundaryEdges"]


(* ::Subsubsection::Closed:: *)
(*Cell Graph*)


Options[ShowCellGraph] = {
	CellVertexStyle -> Directive[Black, EdgeForm[None]],
	ShowVertexLabels -> True,
	ShowIntraCellEdges -> True,
	ShowInterCellEdges -> True,
	IntraCellEdgeStyle -> Blue,
	InterCellEdgeStyle -> Red,
	ShowEdgeTranslations -> False
};
ShowCellGraph[cgraph_, opts:OptionsPattern[{ShowCellGraph, Graph}]] := Module[
	{format\[Gamma]},
	format\[Gamma][expr_]:=StringReplace[expr,{RegularExpression["g(\\d+)(\\^(\\-?\\d+))?"]->ToString[StringForm[\!\(\*
TagBox[
StyleBox["\"\<\\!\\(\\*SuperscriptBox[SubscriptBox[\\(\\[Gamma]\\), \\(`1`\\)], \\(`2`\\)]\\)\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),"$1","$3"],StandardForm],"*"->""}];
	
	Graph[
		cgraph["Graph"],
		Sequence@@FilterRules[{opts}, Options[Graph]],
		
		(* default options for vertices *)
		VertexSize -> AbsolutePointSize[0.2],
		VertexStyle -> OptionValue[CellVertexStyle],
		VertexLabels -> If[OptionValue[ShowVertexLabels], "Name", None],
		
		(* default options for edges *)
		EdgeLabels -> If[OptionValue[ShowEdgeTranslations],
			Thread[EdgeList@cgraph["Graph"] -> format\[Gamma]/@cgraph["EdgeTranslations"]][[
				Position[cgraph["EdgeTranslations"],#][[1, 1]]&/@Select[cgraph["EdgeTranslations"], #!="1"&]
			]], None],
		EdgeStyle -> Thread[EdgeList@cgraph["Graph"] -> (If[#=="1",
			If[OptionValue[ShowIntraCellEdges], OptionValue[IntraCellEdgeStyle], Opacity[0]],
			If[OptionValue[ShowInterCellEdges], OptionValue[InterCellEdgeStyle], Opacity[0]]]&/@cgraph["EdgeTranslations"])
		],
		EdgeShapeFunction -> GraphElementData[{"FilledArcArrow", "ArrowSize" -> .01}],
		Sequence@@FilterRules[{opts}, Options[Graph]]
	]
]


Options[ShowCellSchwarzTriangles] = {
	ShowLabels -> False,
	TriangleRange -> All
};
ShowCellSchwarzTriangles[cgraph_, opts:OptionsPattern[{ShowCellSchwarzTriangles, Graphics}]] := Show[
	Graphics[{
			GetSchwarzTriangle[cgraph["TriangleGroup"], #]&
			/@cgraph["SchwarzTriangleLabels"][[OptionValue[TriangleRange]]]
		},
		Sequence@@FilterRules[{opts}, Options[Graphics]]
	]
]


Options[ShowCellGraphFlattened] = {	
	CellVertexStyle -> Directive[Black, EdgeForm[None]],
	ShowVertexLabels -> True,
	ShowIntraCellEdges -> True,
	ShowInterCellEdges -> True,
	IntraCellEdgeStyle -> Blue,
	InterCellEdgeStyle -> Red
};
ShowCellGraphFlattened[cgraph_, opts:OptionsPattern[{ShowCellGraphFlattened, Graphics, Graph}]] := Module[
	{intracedges, intercedges},
	intracedges = Select[Transpose[{EdgeList@cgraph["Graph"], cgraph["EdgeTranslations"]}], #[[2]]=="1"&][[;;,1]];
	intercedges = Select[Transpose[{EdgeList@cgraph["Graph"], cgraph["EdgeTranslations"]}], #[[2]]!="1"&][[;;,1]];
	
	Show[
		(* edges *)
		Graphics[{
			{OptionValue[IntraCellEdgeStyle],
				Table[GetEdge[cgraph["TriangleGroup"], ResolveEdge[cgraph,edge]],
					{edge, intracedges}
				]
			},
			{OptionValue[InterCellEdgeStyle],
				Table[GetEdge[cgraph["TriangleGroup"], ResolveEdge[cgraph,edge]],
					{edge, intercedges}
				]
			}
		}, Sequence@@FilterRules[{opts}, Options[Graphics]]], 
		(* vertices *)
		Graph[cgraph["Graph"],
			Sequence@@FilterRules[{opts}, Options[Graph]],
			VertexStyle -> OptionValue[CellVertexStyle],
			VertexSize -> AbsolutePointSize[0.2],
			VertexLabels -> If[OptionValue[ShowVertexLabels], "Name", None],
			EdgeStyle -> Opacity[0]
		],
		Sequence@@FilterRules[{opts}, {ImageSize}]
	]
]


Options[ShowCellBoundary] = {
	CellBoundaryStyle -> {Darker@Red, AbsoluteThickness[2]},
	ShowEdgeIdentification -> False,
	EdgeColorFunction -> (ColorData[97,"ColorList"][[Mod[#,15,1]]]&),
	ShowTranslations -> False,
	ShowTranslationLabels -> True,
	ShowTranslationIndices -> False,
	ShowTranslatedCells -> False,
	TranslatedCellBoundaryStyle -> {Black, AbsoluteThickness[1]}
};
ShowCellBoundary[cgraph_, opts:OptionsPattern[]] := Module[{format\[Gamma], gcellbd},
	(* translation label formatting *)
	format\[Gamma][expr_] := StringReplace[expr, {
		RegularExpression["g(\\d+)(\\^(\\-?\\d+))?"] -> ToString[StringForm[\!\(\*
TagBox[
StyleBox[
RowBox[{"\n", "\t\t\t", "\"\<\\!\\(\\*SuperscriptBox[SubscriptBox[\\(\\[Gamma]\\), \\(`1`\\)], \\(`2`\\)]\\)\>\""}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),
			"$1", "$3"
		], StandardForm],
		"*" -> ""
	}];
	
	(* get cell boundary *)
	gcellbd = GetCellBoundary[cgraph];
	
	(* graphics *)
	Graphics[{
		If[OptionValue[ShowEdgeIdentification],
			{Sequence@@OptionValue[CellBoundaryStyle],
				{OptionValue[EdgeColorFunction][#[[1]]],
					#[[2]],
					With[{
						pt = #[[3]]
					},{
						If[OptionValue[ShowTranslations], pt, {}],
						If[OptionValue[ShowTranslationLabels], Text[format\[Gamma]@#[[4]],pt[[1]]], {}],
						If[OptionValue[ShowTranslationIndices], Text[#[[1]],pt[[1]]],{}] 
					}]}&/@gcellbd
			},
			{{Sequence@@OptionValue[CellBoundaryStyle],
				LToGraphics[LLine[{
					GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[1,1]],#[[1]]],
					GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[2,1]],#[[2]]]
				}]&/@cgraph["BoundaryEdges"], Model -> PoincareDisk]
			}, With[{
				pt = LToGraphics[GetSitePosition[cgraph["TriangleGroup"], 3,
					ResolveTranslation[#[[6]], cgraph["TranslationGenerators"]]],
					Model->PoincareDisk
				]},{
					If[OptionValue[ShowTranslations], pt, {}],
					If[OptionValue[ShowTranslationLabels], Text[format\[Gamma]@#[[6]],pt[[1]]], {}]
				}]&/@cgraph["BoundaryEdges"]
			}
		],
		If[OptionValue[ShowTranslatedCells],
			{Sequence@@OptionValue[TranslatedCellBoundaryStyle],
				LToGraphics[Table[LLine[{
						GetSitePosition[cgraph["TriangleGroup"],
							EdgeList[cgraph["Graph"]][[#[[3]]]][[1,1]],
							#[[1]]<>"*"<>ResolveTranslation[\[Gamma], cgraph["TranslationGenerators"]]
						],
						GetSitePosition[cgraph["TriangleGroup"],
							EdgeList[cgraph["Graph"]][[#[[3]]]][[2,1]],
							#[[2]]<>"*"<>ResolveTranslation[\[Gamma],cgraph["TranslationGenerators"]]
						]
					}]&/@cgraph["BoundaryEdges"],
					{\[Gamma], cgraph["BoundaryEdges"][[;;,6]]}],
					Model -> PoincareDisk
				]
			},
		{}]
	}, Sequence@@FilterRules[{opts}, Options[Graphics]]]
]


Options[VisualizeCellGraph] = {
	Elements -> <|
		ShowCellGraph -> {},
		ShowCellBoundary -> {}
	|>
};
VisualizeCellGraph[cgraph_, opts:OptionsPattern[{VisualizeCellGraph, ShowTriangles, GetTriangleTessellation, Graphics, Rasterize}]] := Module[{
		sel = KeySelect[MemberQ[{
			ShowCellSchwarzTriangles,
			ShowCellBoundary,
			ShowCellGraph,
			ShowCellGraphFlattened
		},#]&]
	},
	
	Show[
		ShowTriangles[cgraph["TriangleGroup"],
			Sequence@@FilterRules[{opts},Join@@(Options/@{
				ShowTriangles,
				GetTriangleTessellation,
				Graphics,
				Rasterize
			})]
		],
		KeyValueMap[#1[cgraph, Sequence@@#2]&, sel[OptionValue[Elements]]],
		Sequence@@FilterRules[{opts}, {ImageSize}]
	]
]


(* ::Subsubsection:: *)
(*Model Cell Graph*)


ResolveDGraphEdge[dgraph_,cgraph_,edge_]:=Module[{\[Gamma],v1,v2},
\[Gamma]=ResolveTranslation[dgraph["EdgeTranslations"][[Position[EdgeList@dgraph["Graph"],edge][[1,1]]]],cgraph["TranslationGenerators"]];

v1=ResolveVertex[dgraph,edge[[1]]];
v2=ResolveVertex[dgraph,edge[[2]]];

{v1,{v2[[1]],v2[[2]]<>"*"<>\[Gamma]}}
]
ResolveDGraphConjEdge[dgraph_,cgraph_,edge_]:=Module[{\[Gamma],v1,v2},
\[Gamma]=ResolveTranslation[dgraph["EdgeTranslations"][[Position[EdgeList@dgraph["Graph"],edge][[1,1]]]],cgraph["TranslationGenerators"]];

v1=ResolveVertex[dgraph,edge[[2]]];
v2=ResolveVertex[dgraph,edge[[1]]];

{v1,{v2[[1]],v2[[2]]<>"*("<>\[Gamma]<>")^-1"}}
]


(* ::Input::Initialization:: *)
ClearAll[VisualizeDGraphEdges]
Options[VisualizeDGraphEdges]={ShowLabels->False,CellVertexStyle->Directive[Black,EdgeForm[None]],CellVertexSize->AbsolutePointSize[0.2],ShowIntraCellEdges->True,ShowInterCellEdges->True,IntraCellEdgeStyle->Directive[Blue,AbsoluteThickness[1]],InterCellEdgeStyle->Directive[Red,AbsoluteThickness[1]],ShowCellBoundary->True,ShowConjugateEdges->False,ConjugateEdgeStyle->Dashed,ShowTriangleTessellation->True,EdgeStyleFunction->({}&),VertexStyleFunction->Automatic};
VisualizeDGraphEdges[dgraph_,cgraph_,opts:OptionsPattern[{VisualizeDGraphEdges,Graph,GCellBoundary}]]:=Module[{intracedges,intercedges,inve},
intracedges=Select[Transpose[{EdgeList@dgraph["Graph"],dgraph["EdgeTranslations"]}],#[[2]]=="1"&][[;;,1]];
intercedges=Select[Transpose[{EdgeList@dgraph["Graph"],dgraph["EdgeTranslations"]}],#[[2]]!="1"&][[;;,1]];

Show[
If[OptionValue[ShowTriangleTessellation],ShowTriangles[dgraph["TriangleGroup"]],{}],
Graphics[{
If[OptionValue[ShowIntraCellEdges],
{OptionValue[IntraCellEdgeStyle],Table[{Sequence[OptionValue[EdgeStyleFunction][edge]],GetEdge[dgraph["TriangleGroup"],ResolveDGraphEdge[dgraph,cgraph,edge]]},{edge,intracedges}]},{}],
If[OptionValue[ShowInterCellEdges],
{OptionValue[InterCellEdgeStyle],Table[{Sequence[OptionValue[EdgeStyleFunction][edge]],GetEdge[dgraph["TriangleGroup"],ResolveDGraphEdge[dgraph,cgraph,edge]]},{edge,intercedges}],
OptionValue[ConjugateEdgeStyle],
If[OptionValue[ShowConjugateEdges],Table[{Sequence[OptionValue[EdgeStyleFunction][edge]],GetEdge[dgraph["TriangleGroup"],ResolveDGraphConjEdge[dgraph,cgraph,edge]]},{edge,intercedges}],{}]
},{}]
}/.Line[a_]:>Line[Re@a]],
Graph[dgraph["Graph"],Sequence@@FilterRules[{opts},Options[Graph]],VertexStyle->If[OptionValue[VertexStyleFunction]===Automatic,OptionValue[CellVertexStyle],#->OptionValue[VertexStyleFunction][#]&/@VertexList@dgraph["Graph"]],VertexSize->OptionValue[CellVertexSize],VertexLabels->If[OptionValue[ShowLabels],"Name",None],EdgeStyle->Opacity[0]],
If[OptionValue[ShowCellBoundary],GCellBoundary[cgraph,Sequence@@FilterRules[{opts},Options[GCellBoundary]],ShowEdgeIdentification->True,CellBoundaryStyle->{AbsoluteThickness[2],Opacity[0.5]}],{}],
ImageSize->500]
]


(* ::Section:: *)
(*Package Footer*)


End[];
EndPackage[];
