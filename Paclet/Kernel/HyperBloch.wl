(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["PatrickMLenggenhager`HyperBloch`"];


HCCellGraph::usage = "HCCellGraph[assoc] represents a cell graph with its properties defined by the Association assoc";
HCModelGraph::usage = "HCModelGraph[assoc] represents a model graph with its properties defined by the Association assoc";
HCSupercellModelGraph::usage = "HCSupercellModelGraph[assoc] represents a supercell model graph with its properties defined by the Association assoc";


ImportCellGraphString::usage = "ImportCellGraphString[\"string\"] imports a cell graph from a string and returns an HCCellGraph";
ImportModelGraphString::usage = "ImportModelGraphString[\"string\"] imports a model graph from a string and returns an HCModelGraph";
ImportSupercellModelGraphString::usage = "ImportSupercellModelGraphString[\"string\"] imports a supercell model graph from a string and returns an HCSupercellModelGraph";

HCExampleData::usage = "HCExampleData[\"name\"] imports and returns the specified HCC/HCM/HCS example file from \"PatrickMLenggenhager/HyperBloch/ExampleData/\".";

ShowTriangles::usage = "ShowTriangles[tg] constructs the Schwarz triangles of the triangle group with signature tg in the Poincar\[EAcute] disk representation";

ResolveVertex;
ResolveEdge;
ResolveTranslation;

GetWyckoffPosition::usage = "GetWyckoffPosition[tg, {w,\"g\"}] returns the Cartesian coordinates of the maximally-symmetric Wyckoff position of type w and label g of the triangle group tg";
GetSchwarzTriangle::usage = "GetSchwarzTriangle[tg, \"g\"] returns a Polygon representing the Schwarz triangle labeled by g of the triangle group tg";
GetVertex::usage = "GetVertex[tg, {w,\"g\"}] returns a Point representing a vertex at the position of the maximally-symmetric Wyckoff position of type w and label g of the triangle group tg";
GetEdge::usage = "GetEdge[tg,{{w1, g1},{w2, g2},...}] returns a Line representing the edge (or succession of edges) specified by vertices {wi, gi} of the triangle group tg";
GetCellGraphVertex::usage = "GetCellGraphVertex[cgraph, vertex] returns a Point representing the position of vertex of the HCCellGraph, HCModelGraph, or HCSupercellModelGraph cgraph in the Poincar\[EAcute] disk";
GetCellGraphEdge::usage = "GetCellGraphEdge[cgraph, edge] returns a Line representing edge of the HCCellGraph, HCModelGraph, or HCSupercellModelGraph cgraph in the Poincar\[EAcute] disk";
GetTranslatedCellGraphEdge::usage = "GetTranslatedCellGraphEdge[cgraph, edge, \"\[Gamma]\"] returns a Line representing edge of the HCCellGraph, HCModelGraph, or HCSupercellModelGraph cgraph in the Poincar\[EAcute] disk, translated by the translation \[Gamma]";
GetCellGraphFace::usage = "GetCellGraphFace[cgraph, face] constructs a polygon and list of arrows representing the face face of the cell, model, or supercell model graph cgraph with head HCCellGraph, HCModelGraph, or HCSupercellModelGraph, respectively, and its boundary in the Poincar\[EAcute] disk, respectively";

GetCellBoundary::usage = "GetCellBoundary[cgraph] constructs the cell boundary of the HCCellGraph cgraph";

ShowCellGraph::usage = "ShowCellGraph[cgraph] shows the cell, model, or supercell model graph cgraph with head HCCellGraph, HCModelGraph, or HCSupercellModelGraph, respectively, in the Poincar\[EAcute] disk";
ShowCellSchwarzTriangles::usage = "ShowCellSchwarzTriangles[cgraph] shows the Schwarz triangles making up the cell underlying the cell, model, or supercell model graph cgraph with head HCCellGraph, HCModelGraph, or HCSupercellModelGraph, respectively";
ShowCellGraphFlattened::usage = "ShowCellGraphFlattened[cgraph] shows a flattened, i.e., not compactified, version of the cell, model, or supercell model graph cgraph with head HCCellGraph, HCModelGraph, or HCSupercellModelGraph, respectively, in the Poincar\[EAcute] disk";
ShowCellBoundary::usage = "ShowCellBoundary[cgraph] shows the boundary and boundary identification of the cell on which the cell graph cgraph with head HCCellGraph is defined";

VisualizeCellGraph::usage = "VisualizeCellGraph[cgraph] visualizes the cell graph cgraph with head HCCellGraph in the Poincar\[EAcute] disk with the Schwarz triangles in the background";
VisualizeModelGraph::usage = "VisualizeModelGraph[mgraph] visualizes the (supercell) model graph mgraph with head HCModelGraph (HCSupercellModelGraph) in the Poincar\[EAcute] disk with the Schwarz triangles in the background";

AbelianBlochHamiltonianExpression::usage = "AbelianBlochHamiltonianExpression[mgraph, norb, onsite, hoppings, k] constructs the Abelian Bloch Hamiltonian \[ScriptCapitalH](k) of the HCModelGraph or HCSupercellModelGraph mgraph with the number of orbitals at each site specified by norb, the onsite term by onsite, and the hopping along an edge by hoppings in terms of momenta k[i]";
AbelianBlochHamiltonian::usage = "AbelianBlochHamiltonian[mgraph, norb, onsite, hoppings] returns the Abelian Bloch Hamiltonian \[ScriptCapitalH](k) of the HCModelGraph or HCSupercellModelGraph mgraph with the number of orbitals at each site specified by norb, the onsite term by onsite, and the hopping along an edge by hoppings as a function k :> \[ScriptCapitalH](k)";

NonReciprocalAbelianBlochHamiltonianExpression::usage = "NonReciprocalAbelianBlochHamiltonianExpression[mgraph, norb, onsite, hoppingsCanonical, hoppingsOpposite, k] constructs the non-reciprocal Abelian Bloch Hamiltonian \[ScriptCapitalH](k) of the HCModelGraph or HCSupercellModelGraph mgraph with the number of orbitals at each site specified by norb, the onsite term by onsite, and the hopping along an edge in the canonical direction by hoppingsCanonical and opposite direction by hoppingsOpposite in terms of momenta k[i]";
NonReciprocalAbelianBlochHamiltonian::usage = "NonReciprocalAbelianBlochHamiltonian[mgraph, norb, onsite, hoppingsCanonical, hoppingsOpposite] returns the non-reciprocal Abelian Bloch Hamiltonian \[ScriptCapitalH](k) of the HCModelGraph or HCSupercellModelGraph mgraph with the number of orbitals at each site specified by norb, the onsite term by onsite, and the hopping along an edge in the canonical direction by hoppingsCanonical and opposite direction by hoppingsOpposite as a function k :> \[ScriptCapitalH](k)";


RasterizeGraphics;
NumberOfGenerations;
DiskCenter;
ColorFill;
ColorBoundary;
LineThickness;

SchwarzTriangleOrientation;

ShowEquivalentEdge;

StartingVertex;

CellBoundaryStyle;
EdgeColorFunction;
ShowEdgeIdentification;
ShowTranslatedCells;
ShowTranslationIndices;
ShowTranslationLabels;
ShowTranslations;
TranslatedCellBoundaryStyle;
TranslationLabelStyle;

CellVertexStyle;
EdgeArrowSize;
EdgeFilter;
InterCellEdgeStyle;
IntraCellEdgeStyle;
ShowEdgeTranslations;
ShowInterCellEdges;
ShowIntraCellEdges;
ShowVertexLabels;

CellEdgeStyle;

ShowTriangleLabels;
TriangleLabelStyle;
TriangelRange;
TriangleStyle;

Elements;
CellGraph;

PCModel;
ReturnSparseArray;
CompileFunction;
Parameters;
PBCCluster;


Begin["`Private`"];


(* ::Section:: *)
(*Initialization*)


(* print banner *)
Print["HyperBloch - Version 0.9.1\nMain author: Patrick M. Lenggenhager\n\nThis package loads the following dependencies:\n\t- L2Primitives by Srdjan Vukmirovic\n\t- NCAlgebra by J. William Helton and Mauricio de Oliveira"];


Needs["PatrickMLenggenhager`HyperBloch`L2Primitives`"];


(* check NCAlgebra *)
HyperBloch::NCAlgebra = "`1`";
If[Not[Or@@((Count[ToExpression/@StringCases[#["Version"],
		RegularExpression["(\\d)\\.\\d\\.\\d"]:>"$1"],m_/;m>=6]>0)&/@
		PacletFind["NCAlgebra"])
	],
	Message[HyperBloch::NCAlgebra,
		"No paclet version of NCAlgebra with version >= 6 was found. Version 6 or later is required. Custom installations should work but might lead to problems."
	];
];


(* load NCAlgebra *)
Quiet[Get["NCAlgebra`"], {NCAlgebra`NCAlgebra::SmallCapSymbolsNonCommutative}];
SetOptions[inv, Distribute->True];
SetCommutative[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z];
ClearAll[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z];


(* ::Section:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*Type Definitions*)


HCCellGraph[cgraph_][key_] := cgraph[key]
HCModelGraph[mgraph_][key_] := mgraph[key]
HCSupercellModelGraph[scmgraph_][key_] := scmgraph[key]


(* ::Subsection::Closed:: *)
(*Helper Functions*)


GetFullGraph[graph_] := Graph[
	VertexList@graph,
	EdgeList@DirectedGraph[EdgeList@graph /. DirectedEdge -> UndirectedEdge],
	Evaluate@Sequence@@AbsoluteOptions[graph, VertexCoordinates] 
]


GetUndirectedGraph[graph_] := Graph[
	VertexList@graph,
	EdgeList@Graph[EdgeList@graph /. DirectedEdge -> UndirectedEdge],
	Evaluate@Sequence@@AbsoluteOptions[graph, VertexCoordinates]
]


CyclicallyPermuteList[list_, n_] := Permute[list, PermutationPower[Cycles[{Range[Length[list]]}], n]]


CyclicallyPermuteFaceEdges[face_, n_]:=Graph[VertexList[face], CyclicallyPermuteList[EdgeList[face], n], Sequence@@AbsoluteOptions[face]]


(* ::Subsection::Closed:: *)
(*Group Elements and Vertex Positions*)


InterpretGroupElementString[str_, rules_] := NCExpand@ToExpression@StringReplace[
	str,
	Join[{RegularExpression["\\^([\\-0-9]+)"]->"^($1)","*"->"**"}, (#[[1]] -> ToString[#[[2]]])&/@ rules]
];


Options[GetSitePosition] = {SchwarzTriangleOrientation -> "Default", DiskCenter -> 3};
GetSitePosition[tg_, w_, expr_, OptionsPattern[]] := Module[
	{p, q, r, P, Q, R, ops, triangle, op, i, rules},
	(* triangle group signature *)
	{r, q, p} = Sort[tg];
	
	(* position of sites of fundamental triangle *)
	{R, Q, P} = Switch[OptionValue[SchwarzTriangleOrientation],
		"Default", Block[{pp, qq, rr, pts}, 
			{rr, qq, pp} = CyclicallyPermuteList[{r, q, p}, Mod[-OptionValue[DiskCenter], 3]];
			pts = {
				PDPoint[\[Sqrt]((Cos[\[Pi] * (1/pp+1/rr)]+Cos[\[Pi]/qq])/(Cos[\[Pi] * (1/pp-1/rr)]+Cos[\[Pi]/qq])){Cos[Pi/pp], Sin[Pi/pp]}],
				PDPoint[\[Sqrt]((Cos[\[Pi] * (1/pp+1/qq)]+Cos[\[Pi]/rr])/(Cos[\[Pi] * (1/pp-1/qq)]+Cos[\[Pi]/rr])){1,0}],
				PDPoint[{0,0}]
			};
			CyclicallyPermuteList[pts, -Mod[-OptionValue[DiskCenter], 3]]
		]
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
		})[{R, Q, P}[[w]]]
	]
]


(* ::Subsection::Closed:: *)
(*Import of Cell Graphs*)


ImportCellGraphString[str_]:=Module[{
		version, tg, specs, rels, center, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vlbls, vcoords,
		vertices, edges, etransls, facesstr, faces, faceedges, boundary
	},
	If[StringStartsQ[str, "HyperCells"],
		{version, tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vertices, edges, etransls, facesstr, boundary} =
			StringSplit[StringReplace[str,{"["->"{","]"->"}"}], "\n"];
		version = StringReplace[version, RegularExpression["HyperCells HCC version ([0-9.]+)"] -> "$1"];,
		{tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, vertices, edges, etransls, facesstr, boundary} =
			StringSplit[StringReplace[str,{"["->"{","]"->"}"}], "\n"];
		version = "";
	];

	(* info *)
	tg = ToExpression@tg; (* triangle group signature *)
	{rels, center} = StringTrim[#, {"{", "}"}]&/@ StringSplit[specs, "},"];
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
					DirectedEdge[edges[[e[[1]], 2]], edges[[e[[1]], 1]], -edges[[e[[1]], 3]]]
				],
				{e, face}
			]
		],
		{c, 1, 3}, {face, ToExpression[facesstr][[c]]}
	];
	faceedges = Map[edges[[#[[1]]]]&,ToExpression[facesstr][[center]], {2}];
	
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
		LToGraphics[GetSitePosition[tg, vertices[[i,1]], vlbls[[i]], DiskCenter -> center], Model->PoincareDisk][[1]],
		{i, Length[vertices]}
	];
	
	HCCellGraph[<|
		"TriangleGroup" -> tg,
		"CellCenter" -> center,
		"Genus" -> Length[\[CapitalGamma]gens]/2,
		"Graph"-> Graph[vertices, edges, VertexCoordinates -> vcoords],
		"VertexLabels" -> vlbls,
		"SchwarzTriangleLabels" -> TD\[CapitalGamma],
		"EdgeTranslations" -> etransls,
		"BoundaryEdges" -> boundary,
		"TranslationGenerators" -> \[CapitalGamma]gens,
		"Faces" -> faces[[center]],
		"FaceEdges" -> faceedges,
		"AllFaces" -> faces
	|>]
]


(* ::Subsection::Closed:: *)
(*Import of Model Graphs*)


ImportModelGraphString[str_]:=Module[{
		version, tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, model, vertices, edges, etransls, facesstr,
		rels, center, vlbls, vcoords, graph, faces, faceedges
	},
	If[StringStartsQ[str, "HyperCells"],
		{version, tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, model, vertices, edges, etransls, facesstr} =
			StringSplit[StringReplace[str, {"["->"{", "]"->"}"}], "\n"];
		version = StringReplace[version, RegularExpression["HyperCells HCM version ([0-9.]+)"] -> "$1"];,
		{tg, specs, \[CapitalGamma]gens, TD\[CapitalGamma], TGGw, model, vertices, edges, etransls, facesstr} =
			StringSplit[StringReplace[str, {"["->"{", "]"->"}"}], "\n"];
		version = "";
	];
	
	
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
					DirectedEdge[edges[[e[[1]], 2]], edges[[e[[1]], 1]], -edges[[e[[1]], 3]]]
				],
				{e, face}
			]
		],
		{face, ToExpression[facesstr]}
	];
	faceedges = Map[edges[[#[[1]]]]&, ToExpression[facesstr], {2}];
	
	(* vertex labels and coordinates *)
	vlbls = (StringSplit[StringTrim[#, {"{ ", " }"}], ", "]&/@StringSplit[StringTrim[TGGw, {"{ ", " }"}], " }, { "])[[#[[1]], #[[2]]]]&/@vertices;	
	vcoords = Table[
		LToGraphics[GetSitePosition[tg, vertices[[i,1]], vlbls[[i]], DiskCenter -> center], Model->PoincareDisk][[1]],
		{i, Length[vertices]}
	];
	
	(* edge translations *)
	etransls = StringTrim[StringSplit[StringTrim[etransls,{"{", "}"}], ","], " "];
	
	(* graph *)
	graph = Graph[vertices, edges, VertexCoordinates -> vcoords];
	
	HCModelGraph[<|
		"TriangleGroup" -> tg,
		"CellCenter" -> center,
		"Genus" -> Length[\[CapitalGamma]gens]/2,
		"Graph" -> graph,
		"UndirectedGraph" -> GetUndirectedGraph@graph,
		"FullGraph" -> GetFullGraph@graph,
		"VertexLabels" -> vlbls,
		"SchwarzTriangleLabels" -> TD\[CapitalGamma],
		"EdgeTranslations" -> etransls,
		"TranslationGenerators" -> \[CapitalGamma]gens,
		"Faces" -> faces,
		"FaceEdges" -> faceedges
	|>]
]


(* ::Subsection::Closed:: *)
(*Import of Supercell Model Graphs*)


ImportSupercellModelGraphString[str_]:=Module[{
		version, tg, specs, \[CapitalGamma]0gens, TD\[CapitalGamma]0, TG0Gw, \[CapitalGamma]gens, \[CapitalGamma]\[CapitalGamma]0, T\[CapitalGamma]0\[CapitalGamma], TD\[CapitalGamma], TGGw,
		model, vertices, vertexpos, edges, etransls, facesstr,
		rels0, rels, center, vlbls, vcoords, graph, faces, faceedges
	},
	
	If[StringStartsQ[str, "HyperCells"],
		{version, tg, specs, \[CapitalGamma]0gens, TD\[CapitalGamma]0, TG0Gw, \[CapitalGamma]gens, \[CapitalGamma]\[CapitalGamma]0, T\[CapitalGamma]0\[CapitalGamma], TD\[CapitalGamma], TGGw, model, vertices, vertexpos, edges, etransls, facesstr} =
			StringSplit[StringReplace[str, {"["->"{", "]"->"}"}], "\n"];
		version = StringReplace[version, RegularExpression["HyperCells HCS version ([0-9.]+)"] -> "$1"];,
		{tg, specs, \[CapitalGamma]0gens, TD\[CapitalGamma]0, TG0Gw, \[CapitalGamma]gens, \[CapitalGamma]\[CapitalGamma]0, T\[CapitalGamma]0\[CapitalGamma], TD\[CapitalGamma], TGGw, model, vertices, vertexpos, edges, etransls, facesstr} =
			StringSplit[StringReplace[str, {"["->"{", "]"->"}"}], "\n"];
		version = "";
	];
	
	(* info *)
	tg = ToExpression@tg; (* triangle group signature *)
	{rels0, rels, center} = StringTrim[#, {"{", "}"}]&/@ StringSplit[specs, "},"];
	rels0 = StringSplit[rels0, ", "]; (* primitive cell relators *)
	rels = StringSplit[rels, ", "]; (* supercell relators *)
	center = ToExpression@center; (* cell center *)
	
	(* algebra *)
	\[CapitalGamma]0gens = "(" <> # <> ")" &/@(AssociationThread@@(StringTrim[StringSplit[StringTrim[#, {"{", "}"}], ","], " "]&/@StringSplit[\[CapitalGamma]0gens, " -> "]));
	TD\[CapitalGamma]0 = StringTrim[StringSplit[StringTrim[TD\[CapitalGamma]0, {"{", "}"}], ","], " "];
	\[CapitalGamma]gens = "(" <> # <> ")" &/@(AssociationThread@@(StringTrim[StringSplit[StringTrim[#, {"{", "}"}], ","], " "]&/@StringSplit[\[CapitalGamma]gens, " -> "]));
	\[CapitalGamma]\[CapitalGamma]0 = "(" <> # <> ")" &/@(AssociationThread@@(StringTrim[StringSplit[StringTrim[#, {"{", "}"}], ","], " "]&/@StringSplit[\[CapitalGamma]\[CapitalGamma]0, " -> "]));
	T\[CapitalGamma]0\[CapitalGamma] = StringTrim[StringSplit[StringTrim[T\[CapitalGamma]0\[CapitalGamma], {"{", "}"}], ","], " "];
	TD\[CapitalGamma] = StringTrim[StringSplit[StringTrim[TD\[CapitalGamma], {"{", "}"}], ","], " "];
	
	(* graph *)
	vertices = ToExpression@vertices;
	vertexpos = StringTrim[StringSplit[StringTrim[vertexpos, {"{", "}"}], ","], " "];
	edges = DirectedEdge[vertices[[#1]], vertices[[#2]], #3]&@@@ToExpression[edges];
	faces = Table[
		Graph[
			Table[
				If[e[[2]] == 1,
					edges[[e[[1]]]],
					DirectedEdge[edges[[e[[1]], 2]], edges[[e[[1]], 1]], -edges[[e[[1]], 3]]]
				],
				{e, face}
			]
		],
		{face, ToExpression[facesstr]}
	];
	faceedges = Map[edges[[#[[1]]]]&, ToExpression[facesstr], {2}];
	
	(* vertex labels and coordinates *)
	vlbls = vertexpos;(*(StringSplit[StringTrim[#, {"{ ", " }"}], ", "]&/@StringSplit[StringTrim[TGGw, {"{ ", " }"}], " }, { "])[[#[[1]], #[[2]]]]&/@vertices;	*)
	vcoords = Table[
		LToGraphics[GetSitePosition[tg, vertices[[i,1]], vlbls[[i]], DiskCenter -> center], Model->PoincareDisk][[1]],
		{i, Length[vertices]}
	];
	
	(* edge translations *)
	etransls = StringTrim[StringSplit[StringTrim[etransls,{"{", "}"}], ","], " "];
	
	(* graph *)
	graph = Graph[vertices, edges, VertexCoordinates -> vcoords];
	
	HCSupercellModelGraph[<|
		"TriangleGroup" -> tg,
		"CellCenter" -> center,
		"PCGenus" -> Length[\[CapitalGamma]0gens]/2,
		"Genus" -> Length[\[CapitalGamma]gens]/2,
		"Graph" -> graph,
		"UndirectedGraph"->GetUndirectedGraph@graph,
		"FullGraph" -> GetFullGraph@graph,
		"VertexLabels" -> vlbls,
		"SchwarzTriangleLabels" -> TD\[CapitalGamma],
		"EdgeTranslations" -> etransls,
		"PCTranslationGenerators" -> \[CapitalGamma]0gens,
		"TranslationGenerators" -> \[CapitalGamma]gens,
		"TranslationGroupEmbedding" -> \[CapitalGamma]\[CapitalGamma]0,
		"InternalSupercellTranslations" -> T\[CapitalGamma]0\[CapitalGamma],
		"Faces" -> faces,
		"FaceEdges" -> faceedges
	|>]
]


(* ::Subsection::Closed:: *)
(*Example Data*)


HCExampleData[filename_] := Module[{
	content = Import["PatrickMLenggenhager/HyperBloch/ExampleData/"<>filename]
},
	Switch[StringSplit[filename, "."][[-1]],
		"hcc", ImportCellGraphString[content],
		"hcm", ImportModelGraphString[content],
		"hcs", ImportSupercellModelGraphString[content]
	]
]


(* ::Subsection:: *)
(*Graphical Visualization*)


(* ::Subsubsection::Closed:: *)
(*Triangle Tessellations*)


(* Tomas' code *)
ClearAll[GetTriangleTessellation];
Options[GetTriangleTessellation]={ColorBoundary->RGBColor[{0,0,0,0}],ColorFill->LightGray,LineThickness->0};
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
Q=PDPoint[Sqrt[(Cos[\[Pi] * (1/p+1/q)]+Cos[\[Pi]/r])/(Cos[\[Pi] * (1/p-1/q)]+Cos[\[Pi]/r])]*{1,0}];
R=PDPoint[Sqrt[(Cos[\[Pi] * (1/p+1/r)]+Cos[\[Pi]/q])/(Cos[\[Pi] * (1/p-1/r)]+Cos[\[Pi]/q])]*{Cos[Pi/p],Sin[Pi/p]}];



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
	NumberOfGenerations -> 2,
	DiskCenter -> 3
};
ShowTriangles[tg_, opts:OptionsPattern[{ShowTriangles, Graphics, Rasterize, GetTriangleTessellation}]] := ShowTriangles[tg, opts] = Graphics[
	{
		If[OptionValue[RasterizeGraphics],
			Inset[Rasterize[
					Show[
						GetTriangleTessellation[CyclicallyPermuteList[Sort@tg, Mod[-OptionValue[DiskCenter], 3]], OptionValue[NumberOfGenerations],
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
			GetTriangleTessellation[CyclicallyPermuteList[Sort@tg, Mod[-OptionValue[DiskCenter], 3]], OptionValue[NumberOfGenerations],
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


GetWyckoffPosition[tg_, {w_, g_}, opts:OptionsPattern[{GetSitePosition}]] := LToGraphics[
	GetSitePosition[tg, w, g, Sequence@@FilterRules[{opts}, Options[GetSitePosition]]],
	Model -> PoincareDisk
][[1]]
GetSchwarzTriangle[tg_, g_, opts:OptionsPattern[{GetSitePosition}]] := LToGraphics[
	LPolygon[GetSitePosition[tg, #, g, Sequence@@FilterRules[{opts}, Options[GetSitePosition]]]&/@{1, 2, 3}],
	Model -> PoincareDisk
]
GetVertex[tg_, v_, opts:OptionsPattern[{GetSitePosition}]] := LToGraphics[
	GetSitePosition[tg, v[[1]], v[[2]], Sequence@@FilterRules[{opts}, Options[GetSitePosition]]],
	Model -> PoincareDisk
]
GetEdge[tg_, e_, opts:OptionsPattern[{GetSitePosition}]] := LToGraphics[
	LLine[GetSitePosition[tg, #[[1]], #[[2]], Sequence@@FilterRules[{opts}, Options[GetSitePosition]]]&/@e],
	Model -> PoincareDisk
]


ResolveVertex[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, vertex_] := {
	vertex[[1]],
	cgraph["VertexLabels"][[Position[VertexList@cgraph["Graph"], vertex][[1, 1]]]]
}


ResolveTranslation[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, transl_] := StringReplace[transl,
	RegularExpression["g(\\d+)"] :> cgraph["TranslationGenerators"]["g$1"]
]


ResolveEdge[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, edge_] := Module[{\[Gamma], v1, v2},
	If[MemberQ[EdgeList@cgraph["Graph"], edge],
		(* edge with default orientation *)
		\[Gamma] = ResolveTranslation[cgraph,
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge][[1, 1]]]]
		];
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];,
		
		(* edge with inverted orientation *)
		\[Gamma] = "(" <> ResolveTranslation[cgraph,
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge[[{2, 1, 3}]]][[1, 1]]]]
		] <> ")^(-1)";
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];
	];
	
	{v1, {v2[[1]], v2[[2]]<>"*"<>\[Gamma]}}
]


ResolveEquivalentEdge[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, edge_] := Module[{\[Gamma]inv, v1, v2},
	If[MemberQ[EdgeList@cgraph["Graph"], edge],
		(* edge with default orientation *)
		\[Gamma]inv = "(" <> ResolveTranslation[cgraph,
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge][[1, 1]]]]
		] <> ")^(-1)";
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];,
		
		(* edge with inverted orientation *)
		\[Gamma]inv = ResolveTranslation[cgraph,
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge[[{2, 1, 3}]]][[1, 1]]]]
		];
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];
	];
	
	{{v1[[1]], v1[[2]] <> "*" <> \[Gamma]inv}, v2}
]


GetCellGraphVertex[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, vertex_] := GetVertex[
	cgraph["TriangleGroup"],
	ResolveVertex[cgraph, vertex],
	DiskCenter -> cgraph["CellCenter"]
]
Options[GetCellGraphEdge] = { ShowEquivalentEdge -> False };
GetCellGraphEdge[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, edge_, OptionsPattern[]] :=
If[OptionValue[ShowEquivalentEdge],
	GetEdge[
		cgraph["TriangleGroup"],
		ResolveEquivalentEdge[cgraph, edge],
		DiskCenter -> cgraph["CellCenter"]
	],
	GetEdge[
		cgraph["TriangleGroup"],
		ResolveEdge[cgraph, edge],
		DiskCenter -> cgraph["CellCenter"]
	]
]


ResolveTranslatedEdge[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, edge_, \[Gamma]0_] := Module[{\[Gamma], v1, v2},
	If[NumericQ@edge[[3]] && edge[[3]] > 0 || ListQ@edge[[3]] && edge[[3, 1]] > 0,
		(* edge with default orientation *)
		\[Gamma] = ResolveTranslation[cgraph,
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"], edge][[1, 1]]]]
		];
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];,
		(* edge with inverted orientation *)
		\[Gamma] = "(" <> ResolveTranslation[cgraph,
			cgraph["EdgeTranslations"][[Position[EdgeList@cgraph["Graph"],
				DirectedEdge[edge[[2]], edge[[1]], -edge[[3]]]	
			][[1, 1]]]]
		] <> ")^(-1)";
	
		v1 = ResolveVertex[cgraph, edge[[1]]];
		v2 = ResolveVertex[cgraph, edge[[2]]];
	];
	
	{{{v1[[1]], v1[[2]]<>"*"<>\[Gamma]0}, {v2[[1]], v2[[2]]<>"*"<>\[Gamma]<>"*"<>\[Gamma]0}}, \[Gamma]<>"*"<>\[Gamma]0}
]
GetTranslatedCellGraphEdge[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, edge_, \[Gamma]_] := GetEdge[
	cgraph["TriangleGroup"],
	ResolveTranslatedEdge[cgraph, edge, ResolveTranslation[cgraph, \[Gamma]]][[1]],
	DiskCenter -> cgraph["CellCenter"]
]


ResolveFace[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, face_] := Module[{\[Gamma] = "1", re},
	Table[
		re = ResolveTranslatedEdge[cgraph, e, \[Gamma]];
		\[Gamma] = re[[2]];
		re[[1]],
		{e, EdgeList[face]}
	]
]
Options[GetCellGraphFace] = {
	StartingVertex -> "LowestSchwarzTriangleIndex"
};
GetCellGraphFace[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph, face_, opts:OptionsPattern[]] := GetCellGraphFace[cgraph, face, opts] = Module[{
		pos, index
	},
	index = If[ListQ@OptionValue[StartingVertex],
		Switch[OptionValue[StartingVertex][[1]],
		"Index", OptionValue[StartingVertex][[2]],
		"VertexCriterion", FirstPosition[EdgeList[face][[;;,1]], Select[EdgeList[face][[;;,1]], OptionValue[StartingVertex][[2]]][[1]]][[1]],
		"EdgeCriterion", FirstPosition[EdgeList[face], Select[EdgeList[face], OptionValue[StartingVertex][[2]]][[1]]][[1]],
		"VertexPattern", FirstPosition[EdgeList[face][[;;,1]], OptionValue[StartingVertex][[2]]][[1]],
		"EdgePattern", FirstPosition[EdgeList[face], OptionValue[StartingVertex][[2]]][[1]],
		_, 1
		],
		Switch[OptionValue[StartingVertex],
		"LowestSchwarzTriangleIndex", FirstPosition[EdgeList[face][[;;,1]], Min[EdgeList[face][[;;, 1, 2]]]][[1]],
		_, 1
		]
	];
	pos = GetSitePosition[cgraph["TriangleGroup"], #[[1,1]], #[[1,2]], DiskCenter -> cgraph["CellCenter"]]&/@
		ResolveFace[cgraph, CyclicallyPermuteFaceEdges[face, -(index-1)]];
	{
		LToGraphics[LPolygon[pos], Model -> PoincareDisk],
		Arrow@LToGraphics[LLine[{pos[[#]], pos[[Mod[# + 1, Length@pos, 1]]]}], Model -> PoincareDisk]&/@Range[1, Length@pos ]
	}
]


(* ::Subsubsection::Closed:: *)
(*Cell Boundary*)


GetCellBoundary[cgraph_HCCellGraph] := GetCellBoundary[cgraph] = {
	#[[4]],
	LToGraphics[LLine[{
		GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[1, 1]], #[[1]], DiskCenter -> cgraph["CellCenter"]],
		GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[2, 1]], #[[2]], DiskCenter -> cgraph["CellCenter"]]
	}], Model -> PoincareDisk],
	LToGraphics[GetSitePosition[cgraph["TriangleGroup"], 3,
		ResolveTranslation[cgraph, #[[6]]],
		DiskCenter -> cgraph["CellCenter"]
	], Model -> PoincareDisk],
	#[[6]]
}&/@cgraph["BoundaryEdges"]


(* ::Subsubsection::Closed:: *)
(*Cell Graph*)


Options[ShowCellGraph] = {
	CellVertexStyle -> Directive[Black, AbsolutePointSize[5]],
	ShowVertexLabels -> True,
	EdgeArrowSize -> .015,
	ShowIntraCellEdges -> True,
	ShowInterCellEdges -> True,
	IntraCellEdgeStyle -> Blue,
	InterCellEdgeStyle -> Red,
	ShowEdgeTranslations -> False,
	EdgeFilter -> (True&)
};
ShowCellGraph[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph,
	opts:OptionsPattern[{ShowCellGraph, Graph}]] := ShowCellGraph[cgraph, opts] = Module[
	{format\[Gamma]},
	format\[Gamma][expr_]:=StringReplace[expr,{RegularExpression["g(\\d+)(\\^(\\-?\\d+))?"]->ToString[StringForm[\!\(\*
TagBox[
StyleBox["\"\<\\!\\(\\*SuperscriptBox[SubscriptBox[\\(\\[Gamma]\\), \\(`1`\\)], \\(`2`\\)]\\)\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),"$1","$3"],StandardForm],"*"->""}];
	
	Show[
		Graph[
			cgraph["Graph"],
			Sequence@@FilterRules[{opts}, Options[Graph]],
			
			(* default options for vertices *)
			VertexSize -> 0,
			VertexStyle -> Directive[EdgeForm[Opacity[0]],FaceForm[Opacity[0]]],
			VertexLabels -> If[OptionValue[ShowVertexLabels], "Name", None],
			
			(* default options for edges *)
			EdgeLabels -> If[OptionValue[ShowEdgeTranslations],
				Select[Thread[EdgeList@cgraph["Graph"] -> format\[Gamma]/@cgraph["EdgeTranslations"]][[
					Position[cgraph["EdgeTranslations"],#][[1, 1]]&/@Select[cgraph["EdgeTranslations"], #!="1"&]
				]], OptionValue[EdgeFilter][#[[1]]]&], None],
			EdgeStyle -> Join[
				Select[Thread[EdgeList@cgraph["Graph"] -> (If[#=="1",
				If[OptionValue[ShowIntraCellEdges], OptionValue[IntraCellEdgeStyle], Opacity[0]],
				If[OptionValue[ShowInterCellEdges], OptionValue[InterCellEdgeStyle], Opacity[0]]]&/@cgraph["EdgeTranslations"])
			], OptionValue[EdgeFilter][#[[1]]]&],
				Select[Thread[EdgeList@cgraph["Graph"] -> Opacity[0]], Not@OptionValue[EdgeFilter][#[[1]]]&]
			],
			EdgeShapeFunction -> GraphElementData[{"FilledArcArrow", "ArrowSize" -> OptionValue[EdgeArrowSize]}],
			Sequence@@FilterRules[{opts}, Options[Graph]]
		],
		Graphics[{OptionValue[CellVertexStyle],
			Point/@(VertexCoordinates/.AbsoluteOptions[cgraph["Graph"]])
		}]
	]
]


Options[ShowCellSchwarzTriangles] = {
	TriangleStyle -> Directive[FaceForm[Black], EdgeForm[None]],
	ShowTriangleLabels -> False,
	TriangleLabelStyle -> White,
	TriangleRange -> All
};
ShowCellSchwarzTriangles[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph,
	opts:OptionsPattern[{ShowCellSchwarzTriangles, Graphics}]] := ShowCellSchwarzTriangles[cgraph, opts] = Module[
{format},
	format[expr_] := StringReplace[expr,{RegularExpression["([xyz\\)])(\\^(\\-?\\d+))?"]
		-> ToString[StringForm[\!\(\*
TagBox[
StyleBox["\"\<\\!\\(\\*SuperscriptBox[\\(`1`\\), \\(`2`\\)]\\)\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\), "$1", "$3"],StandardForm], "*"->""}];
	Show[
		Graphics[{
			{OptionValue[TriangleStyle],
				With[{st = GetSchwarzTriangle[cgraph["TriangleGroup"], #, DiskCenter -> cgraph["CellCenter"]]},
					{st, If[OptionValue[ShowTriangleLabels],
						Text[Style[format@#, OptionValue[TriangleLabelStyle]], RegionCentroid@st], {}]}
				]&/@cgraph["SchwarzTriangleLabels"][[OptionValue[TriangleRange]]]}
			},
			Sequence@@FilterRules[{opts}, Options[Graphics]]
		]
	]
]


Options[ShowCellGraphFlattened] = {	
	CellVertexStyle -> Directive[Black, AbsolutePointSize[5]],
	ShowVertexLabels -> True,
	CellEdgeStyle -> Arrowheads[{{Small,0.5}}],
	ShowIntraCellEdges -> True,
	ShowInterCellEdges -> True,
	IntraCellEdgeStyle -> Blue,
	InterCellEdgeStyle -> Red,
	EdgeFilter -> (True&)
};
ShowCellGraphFlattened[cgraph_HCCellGraph|cgraph_HCModelGraph|cgraph_HCSupercellModelGraph,
	opts:OptionsPattern[{ShowCellGraphFlattened, Graphics, Graph}]] := ShowCellGraphFlattened[cgraph, opts] = Module[
	{intracedges, intercedges},
	intracedges = Select[Transpose[{EdgeList@cgraph["Graph"], cgraph["EdgeTranslations"]}], OptionValue[EdgeFilter][#[[1]]]&&#[[2]]=="1"&][[;;,1]];
	intercedges = Select[Transpose[{EdgeList@cgraph["Graph"], cgraph["EdgeTranslations"]}], OptionValue[EdgeFilter][#[[1]]]&&#[[2]]!="1"&][[;;,1]];
	
	Show[
		(* edges *)
		Graphics[{OptionValue[CellEdgeStyle],
			If[OptionValue[ShowIntraCellEdges],
				{OptionValue[IntraCellEdgeStyle],
					Table[Arrow@GetEdge[cgraph["TriangleGroup"], ResolveEdge[cgraph, edge], DiskCenter -> cgraph["CellCenter"]],
						{edge, intracedges}
					]
				}, {}],
			If[OptionValue[ShowInterCellEdges],
				{OptionValue[InterCellEdgeStyle],
					Table[Arrow@GetEdge[cgraph["TriangleGroup"], ResolveEdge[cgraph, edge], DiskCenter -> cgraph["CellCenter"]],
						{edge, intercedges}
					]
				},{}]
		}, Sequence@@FilterRules[{opts}, Options[Graphics]]], 
		(* vertices *)
		Graph[cgraph["Graph"],
			Sequence@@FilterRules[{opts}, Options[Graph]],
			VertexSize -> 0,
			VertexStyle -> Directive[EdgeForm[Opacity[0]],FaceForm[Opacity[0]]],
			VertexLabels -> If[OptionValue[ShowVertexLabels], "Name", None],
			EdgeStyle -> Opacity[0]
		],
		Graphics[{OptionValue[CellVertexStyle],
			Point/@(VertexCoordinates/.AbsoluteOptions[cgraph["Graph"]])
		}],
		Sequence@@FilterRules[{opts}, {ImageSize}]
	]
]


Options[ShowCellBoundary] = {
	CellBoundaryStyle -> Directive[Darker@Red, AbsoluteThickness[2]],
	ShowEdgeIdentification -> False,
	EdgeColorFunction -> (ColorData[97,"ColorList"][[Mod[#,15,1]]]&),
	ShowTranslations -> False,
	ShowTranslationLabels -> True,
	ShowTranslationIndices -> False,
	ShowTranslatedCells -> False,
	TranslatedCellBoundaryStyle -> Directive[Black, AbsoluteThickness[1]],
	TranslationLabelStyle -> {}
};
ShowCellBoundary[cgraph_HCCellGraph, opts:OptionsPattern[]] := ShowCellBoundary[cgraph, opts] = Module[{format\[Gamma], gcellbd},
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
			{OptionValue[CellBoundaryStyle],
				{OptionValue[EdgeColorFunction][#[[1]]],
					#[[2]],
					With[{
						pt = #[[3]]
					},{
						If[OptionValue[ShowTranslations], pt, {}],
						If[OptionValue[ShowTranslationLabels], Text[Style[format\[Gamma]@#[[4]],OptionValue[TranslationLabelStyle]],pt[[1]]], {}],
						If[OptionValue[ShowTranslationIndices], Text[#[[1]],pt[[1]]],{}] 
					}]}&/@gcellbd
			},
			{{OptionValue[CellBoundaryStyle],
				LToGraphics[LLine[{
					GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[1,1]], #[[1]], DiskCenter -> cgraph["CellCenter"]],
					GetSitePosition[cgraph["TriangleGroup"], EdgeList[cgraph["Graph"]][[#[[3]]]][[2,1]], #[[2]], DiskCenter -> cgraph["CellCenter"]]
				}]&/@cgraph["BoundaryEdges"], Model -> PoincareDisk]
			}, With[{
				pt = LToGraphics[GetSitePosition[cgraph["TriangleGroup"], 3,
					ResolveTranslation[cgraph, #[[6]]],
					DiskCenter -> cgraph["CellCenter"]],
					Model->PoincareDisk
				]},{
					If[OptionValue[ShowTranslations], pt, {}],
					If[OptionValue[ShowTranslationLabels], Text[Style[format\[Gamma]@#[[6]],OptionValue[TranslationLabelStyle]],pt[[1]]], {}]
				}]&/@cgraph["BoundaryEdges"]
			}
		],
		If[OptionValue[ShowTranslatedCells],
			{OptionValue[TranslatedCellBoundaryStyle],
				LToGraphics[Table[LLine[{
						GetSitePosition[cgraph["TriangleGroup"],
							EdgeList[cgraph["Graph"]][[#[[3]]]][[1,1]],
							#[[1]]<>"*"<>ResolveTranslation[cgraph, \[Gamma]],
							DiskCenter -> cgraph["CellCenter"]
						],
						GetSitePosition[cgraph["TriangleGroup"],
							EdgeList[cgraph["Graph"]][[#[[3]]]][[2,1]],
							#[[2]]<>"*"<>ResolveTranslation[cgraph, \[Gamma]],
							DiskCenter -> cgraph["CellCenter"]
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
		ShowCellBoundary -> {},
		ShowCellGraph -> {}
	|>
};
VisualizeCellGraph[cgraph_HCCellGraph, opts:OptionsPattern[{VisualizeCellGraph, ShowTriangles, GetTriangleTessellation, Graphics, Rasterize}]] := Module[{
		sel = KeySelect[MemberQ[{
			ShowCellSchwarzTriangles,
			ShowCellBoundary,
			ShowCellGraph,
			ShowCellGraphFlattened
		},#]&]
	},
	
	Show[
		ShowTriangles[cgraph["TriangleGroup"], DiskCenter -> cgraph["CellCenter"],
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


(* ::Subsubsection::Closed:: *)
(*Model Graph*)


Options[VisualizeModelGraph] = {
	Elements -> <|
		ShowCellGraphFlattened -> {}
	|>,
	CellGraph -> None
};
VisualizeModelGraph[mgraph_HCModelGraph|mgraph_HCSupercellModelGraph,
	opts:OptionsPattern[{VisualizeModelGraph, ShowTriangles, GetTriangleTessellation, Graphics, Rasterize}]] :=
Module[{
		sel = KeySelect[MemberQ[{
			ShowCellSchwarzTriangles,
			ShowCellGraph,
			ShowCellGraphFlattened
		}, #]&]
	},
	
	Show[
		ShowTriangles[mgraph["TriangleGroup"], DiskCenter -> mgraph["CellCenter"],
			Sequence@@FilterRules[{opts},Join@@(Options/@{
				ShowTriangles,
				GetTriangleTessellation,
				Graphics,
				Rasterize
			})]
		],
		If[MemberQ[Keys@OptionValue[Elements], ShowCellBoundary] &&
				Not[OptionValue[CellGraph] === None],
			ShowCellBoundary[OptionValue[CellGraph],
				Sequence@@(OptionValue[Elements][ShowCellBoundary])
			],
		{}],
		KeyValueMap[#1[mgraph, Sequence@@#2]&, sel[OptionValue[Elements]]],
		Sequence@@FilterRules[{opts}, {ImageSize}]
	]
]


(* ::Subsection::Closed:: *)
(*Construct Bloch Hamiltonians*)


MakeHermitian[H_] := H + ConjugateTranspose@H


ZeroMatrix[n_] := ConstantArray[0, {n, n}]


Options[AbelianBlochHamiltonianExpression] = {
	PCModel -> None,
	ReturnSparseArray -> False
};
AbelianBlochHamiltonianExpression[model_HCModelGraph|model_HCSupercellModelGraph, norb_, onsite_, hoppings_, k_Symbol,
	OptionsPattern[AbelianBlochHamiltonianExpression]] :=
Module[{dimk, verts, Nverts, edges, htest, Hexpr, PCVertex, PCEdge, H, assumptions},
	(* dimension of Abelian Brillouin zone *)
	dimk = 2*model["Genus"];
	
	(* extract vertices *)
	verts = VertexList@model["Graph"];
	Nverts = Length@verts;
	
	(* extract edges *)
	edges = Transpose[{
		EdgeList@model["Graph"],
		ToExpression@StringReplace[#, RegularExpression["g(\\d+)"] -> "(E^(\[ImaginaryI] " <> ToString[k] <> "[$1]))"]&/@
			model["EdgeTranslations"]
	}];
	
	(* mapping to primitive cell *)
	If[OptionValue[PCModel] === None,
		PCVertex[pcvertex_] := pcvertex;
		PCEdge[pcedge_] := pcedge;,
		PCVertex[scvertex_] := scvertex[[{1, 2}]];
		PCEdge[scedge_] := scedge[[0]][
			VertexList[OptionValue[PCModel]["Graph"]][[scedge[[3, 1]]]],
			VertexList[OptionValue[PCModel]["Graph"]][[scedge[[3, 2]]]],
			scedge[[3,3]]
		];
	];
	
	H = If[norb === 1,
		MakeHermitian@Total[SparseArray[{
			{Position[verts, #1[[2]]][[1, 1]], Position[verts, #1[[1]]][[1, 1]]} -> hoppings[PCEdge@#1]#2
		}, Nverts]&@@@edges] + Total[SparseArray[{
			{Position[verts, #][[1, 1]], Position[verts, #][[1, 1]]} -> onsite[PCVertex@#]
		}, Nverts]&/@verts],
		MakeHermitian@Total[SparseArray`SparseBlockMatrix[Join[
			{{Position[verts, #1[[2]]][[1, 1]], Position[verts, #1[[1]]][[1,1]]} -> hoppings[PCEdge@#1]#2},
			Table[{i, i} -> ZeroMatrix[norb[PCVertex@verts[[i]]]], {i, 1, Length@verts}]
		]]&@@@edges] + SparseArray`SparseBlockMatrix[
			{Position[verts, #][[1, 1]], Position[verts, #][[1, 1]]} -> onsite[PCVertex@#]&/@verts
		]
	];
	assumptions = And@@Table[k[i]\[Element]Reals, {i, 1, dimk}];
	H = Map[Simplify[#, assumptions]&, H, {2}];
	
	If[OptionValue[ReturnSparseArray], SparseArray@H, Normal@H]
]


Options[AbelianBlochHamiltonian] = {
	CompileFunction -> False,
	Parameters -> {},
	PBCCluster -> False
};
AbelianBlochHamiltonian[model_HCModelGraph|model_HCSupercellModelGraph, norb_, onsite_, hoppings_,
	opts:OptionsPattern[{AbelianBlochHamiltonianExpression, AbelianBlochHamiltonian, Compile}]] :=
If[OptionValue[PBCCluster],
	Evaluate[
		AbelianBlochHamiltonianExpression[model, norb, onsite, hoppings, k,
			ReturnSparseArray -> False,
			Evaluate@FilterRules[{opts}, Options[AbelianBlochHamiltonianExpression]]]/.
		Join[Table[k[i] -> 0, {i, 1, 2*model["Genus"]}], OptionValue[Parameters]]
	],
	If[OptionValue[CompileFunction],
		Block[{k},
			Compile[Evaluate@Table[{k[i], _Real}, {i, 1, 2*model["Genus"]}],
				Evaluate[
					AbelianBlochHamiltonianExpression[model, norb, onsite, hoppings, k,
						ReturnSparseArray -> False,
						Evaluate@FilterRules[{opts}, Options[AbelianBlochHamiltonianExpression]]
					]/. OptionValue[Parameters]
				],
				Evaluate@FilterRules[{opts}, Options[Compile]]
			]
		],
		Block[{k},
			Function[Evaluate@Table[Symbol["k" <> ToString@i], {i, 1, 2*model["Genus"]}],
				Evaluate[
					AbelianBlochHamiltonianExpression[model, norb, onsite, hoppings, k,
						ReturnSparseArray -> False,
						Evaluate@FilterRules[{opts}, Options[AbelianBlochHamiltonianExpression]]]/.
					Join[
						Table[k[i] -> Symbol["k" <> ToString@i], {i, 1, 2*model["Genus"]}],
						OptionValue[Parameters]
					]
				]
			]
		]
	]
]


(* ::Subsection::Closed:: *)
(*Construct non-reciprocal Bloch Hamiltonians*)


Options[NonReciprocalAbelianBlochHamiltonianExpression] = {
	PCModel -> None,
	ReturnSparseArray -> False
};
NonReciprocalAbelianBlochHamiltonianExpression[model_HCModelGraph|model_HCSupercellModelGraph, norb_, onsite_, hoppingsCanonical_, hoppingsOpposite_, k_Symbol,
	OptionsPattern[NonReciprocalAbelianBlochHamiltonianExpression]] :=
Module[{dimk, verts, Nverts, edges, htest, Hexpr, PCVertex, PCEdge, H, assumptions},
	(* dimension of Abelian Brillouin zone *)
	dimk = 2*model["Genus"];
	
	(* extract vertices *)
	verts = VertexList@model["Graph"];
	Nverts = Length@verts;
	
	(* extract edges *)
	edges = Transpose[{
		EdgeList@model["Graph"],
		ToExpression@StringReplace[#, RegularExpression["g(\\d+)"] -> "(E^(\[ImaginaryI] " <> ToString[k] <> "[$1]))"]&/@
			model["EdgeTranslations"]
	}];
	
	(* mapping to primitive cell *)
	If[OptionValue[PCModel] === None,
		PCVertex[pcvertex_] := pcvertex;
		PCEdge[pcedge_] := pcedge;,
		PCVertex[scvertex_] := scvertex[[{1, 2}]];
		PCEdge[scedge_] := scedge[[0]][
			VertexList[OptionValue[PCModel]["Graph"]][[scedge[[3, 1]]]],
			VertexList[OptionValue[PCModel]["Graph"]][[scedge[[3, 2]]]],
			scedge[[3,3]]
		];
	];
	
	H = If[norb === 1,
		Total[Normal@SparseArray[{
			{Position[verts, #1[[2]]][[1, 1]], Position[verts, #1[[1]]][[1, 1]]} -> hoppingsCanonical[PCEdge@#1]#2
		}, Nverts]&@@@edges] + ConjugateTranspose@Total[Normal@SparseArray[{
			{Position[verts, #1[[2]]][[1, 1]], Position[verts, #1[[1]]][[1, 1]]} -> Conjugate@hoppingsOpposite[PCEdge@#1]#2
		}, Nverts]&@@@edges] + Total[Normal@SparseArray[{
			{Position[verts, #][[1, 1]], Position[verts, #][[1, 1]]} -> onsite[PCVertex@#]
		}, Nverts]&/@verts],
		Total[Normal@SparseArray`SparseBlockMatrix[Join[
			{{Position[verts, #1[[2]]][[1, 1]], Position[verts, #1[[1]]][[1,1]]} -> hoppingsCanonical[PCEdge@#1]#2},
			Table[{i, i} -> ZeroMatrix[norb[PCVertex@verts[[i]]]], {i, 1, Length@verts}]
		]]&@@@edges] + ConjugateTranspose@Total[Normal@SparseArray`SparseBlockMatrix[Join[
			{{Position[verts, #1[[2]]][[1, 1]], Position[verts, #1[[1]]][[1,1]]} -> Conjugate@hoppingsOpposite[PCEdge@#1]#2},
			Table[{i, i} -> ZeroMatrix[norb[PCVertex@verts[[i]]]], {i, 1, Length@verts}]
		]]&@@@edges] + Normal@SparseArray`SparseBlockMatrix[
			{Position[verts, #][[1, 1]], Position[verts, #][[1, 1]]} -> onsite[PCVertex@#]&/@verts
		]
	];
	assumptions = And@@Table[k[i]\[Element]Reals, {i, 1, dimk}];
	H = Map[Simplify[#, assumptions]&, H, {2}];
	
	If[OptionValue[ReturnSparseArray], SparseArray@H, Normal@H]
]


Options[NonReciprocalAbelianBlochHamiltonian] = {
	CompileFunction -> False,
	Parameters -> {},
	PBCCluster -> False
};
NonReciprocalAbelianBlochHamiltonian[model_HCModelGraph|model_HCSupercellModelGraph, norb_, onsite_,  hoppingsCanonical_, hoppingsOpposite_,
	opts:OptionsPattern[{NonReciprocalAbelianBlochHamiltonianExpression, NonReciprocalAbelianBlochHamiltonian, Compile}]] :=
If[OptionValue[PBCCluster],
	Evaluate[
		NonReciprocalAbelianBlochHamiltonianExpression[model, norb, onsite,  hoppingsCanonical, hoppingsOpposite, k,
			Evaluate@FilterRules[{opts}, Options[NonReciprocalAbelianBlochHamiltonianExpression]]]/.
		Join[Table[k[i] -> 0, {i, 1, 2*model["Genus"]}], OptionValue[Parameters]]
	],
	If[OptionValue[CompileFunction],
		Block[{k},
			Compile[Evaluate@Table[{k[i], _Real}, {i, 1, 2*model["Genus"]}],
				Evaluate[
					NonReciprocalAbelianBlochHamiltonianExpression[model, norb, onsite,  hoppingsCanonical, hoppingsOpposite, k,
						Evaluate@FilterRules[{opts}, Options[NonReciprocalAbelianBlochHamiltonianExpression]]
					]/. OptionValue[Parameters]
				],
				Evaluate@FilterRules[{opts}, Options[Compile]]
			]
		],
		Block[{k},
			Function[Evaluate@Table[Symbol["k" <> ToString@i], {i, 1, 2*model["Genus"]}],
				Evaluate[
					NonReciprocalAbelianBlochHamiltonianExpression[model, norb, onsite,  hoppingsCanonical, hoppingsOpposite, k,
						Evaluate@FilterRules[{opts}, Options[NonReciprocalAbelianBlochHamiltonianExpression]]]/.
					Join[
						Table[k[i] -> Symbol["k" <> ToString@i], {i, 1, 2*model["Genus"]}],
						OptionValue[Parameters]
					]
				]
			]
		]
	]
]


(* ::Section:: *)
(*Package Footer*)


End[];
EndPackage[];
