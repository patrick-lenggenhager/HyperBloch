(* ::Package:: *)

(* : Mathematica Version : 4.0 *)

(* : Name : `L2Primitives`*)

(* : Title : Basic drawing in the various models of the Lobachevskian 
             plane*)

(* : Authors : I. Knezevic, R. Sazdanovic, S. Vukmirovic *)

(* : Copyright : Copyright 2001 - 2002 *)

(* : History : working Version 1.0 at the end of 2001,
               current Version 1.1 at  16.01.2002 *)

(*: Summary : 
      This package provides an internal representation of the basic 
      objects of the Lobachevskian plane, their isometric transforma-
      tions and   presenting (drawing) in three different models
      (KleinDisk, PoincareDisk and   HalfPlane) of the Lobachevskian
       plane. *)

(* : Context : `L2Primitives`*)

(*: Keywords : Lobachevskian plane, Hyperbolic plane *)

(*: Requirements : None.*)

BeginPackage["PatrickMLenggenhager`HyperBloch`L2Primitives`"]; 

   (* globals *)
KleinDisk;
HalfPlane;
PoincareDisk;  
                       (* USAGE MESSAGES *)
                       
Model::usage =
  "Model is an option of LToGraphics specifying the model of the
   Lobachevskian plane in which L2 primitives should be represented.
   Model can be set to: KleinDisk, PoincareDisk, HalfPlane. The 
   default model is KleinDisk in  which LPoint-s are internally 
   represented and isometric transformations are performed."

LPoint::usage =
  "LPoint[{x, y}] internally represents a point in the Lobachevskian
   plane.User is supposed to use PDPoint, KDPoint, HPPoint which will
   return a correct LPoint. Essentially, LPoint has coordinates of 
   a KDPoint."

L2Reflection::usage = 
 "L2Reflection[x_LPoint][y_LPrimitive] returns the result of the 
 reflection of y with respect to x."

L2Reflection::usage =
  "L2Reflection[x_LLine][y_LPrimitive] returns the result of the 
  reflection of y with respect to the line x."

L2Translation::usage = 
"L2Translation[a_LPoint, b_LPoint][c_LPrimitive] returns the result
 of L2 translation (compositon of two L2 reflections with respect to
 lines with common perpendicular ab) of the primitive c."

L2Rotation::usage =
  "L2Rotation[a_LPoint, t_][x_LPrimitive] returns the result of the 
  rotation of x around a by an angle t"

HPPoint::usage =
  "HPPoint[{x, y}] returns LPoint object representing the point (x,y) 
  in the half plane model of L2"
HPPoint::badArg =
  "The second coordinate of HPPoint should be a real number greater
   than zero "

KDPoint::usage =
  "KDPoint[{x, y}] returns LPoint representing the point (x,y) in the 
  KleinDisk model of L2."
KDPoint::badArg =
  "arguments of a KDPoint should be coordinates of a point from the
   unit disk."

PDPoint::usage =
  "PDPoint[{x, y}] returns LPoint representing the point (x,y) in the
   PoincareDisk model of L2."
PDPoint::badArg =
  "arguments of a PDPoint should be coordinates of a point from the 
  unit disk."

LLine::usage =
  "LLine[{x1_LPoint, ..., xn_LPoint}] represents the polygonal line
   connecting the points x1, ..., xn."

LPolygon::usage =
  "LPolygon[{x1_LPoint, ..., xn_LPoint}] represents the polygon dete-
  rmined by the points"

LCircle::usage =
  "LCircle[c_LPoint, r_] represents the circle with L2 centre c and
   L2 radius r" 

LDisk::usage =
  "LDisk[c_LPoint, r_] represents the disk with L2 centre c and L2 
  radius r"

LToGraphics::usage = 
"LToGraphics[x_LPrimitive, opts_] returns Graphics primitive 
corresponding to x in the model specified by an option Model. Option 
Model->model specifies model, option PlotPoints->number specifies
precision (number of segments approximmating the curve) of the 
Graphics object to be created. Models allowed are: PoincareDisk, 
KleinDisk, HalfPlane. LPrimitives are transformed as follows: 
LPoint->Point, LLine->Line, LPolygon->Polygon, LCircle->Circle 
(if model is PoincareDisk, HalfPlane) or Line (if model is KleinDisk),
LDisk->Disk (if model is PoincareDisk, HalfPlane) or Polygon (if
model is KleinDisk). LToGraphics is Listable"

Begin["`Private`"]

Options[LToGraphics] = {Model -> KleinDisk, PlotPoints -> 40};

    (* when to consider line straight in PD or HP model *)
err = 0.001

    (* Euclidean norm of a vector *)
norm[vect_List] := N[Sqrt[vect . vect]]
    (* returns unit vector *)
normalize[vect_List] := If[norm[vect] == 0, 0, vect/norm[vect]]

    (* LPoint (essentially KD point) to projective coordinates *)
LPointToH[LPoint[x_]] := Append[x, 1]
    (* vice versa *)
hToLPoint[x_List] := LPoint[Drop[x, -1]/Last[x]]

    (* Minkowskian scalar product in  proj. coordinates *)
minDot[x_List, y_List] := (x {1, 1, -1}) . y

    (* center of circle given by 3 points *)
CCC[c_List, p1_List, p2_List] := 
  Module[{h},                  
    h = N[normalize[((p2 - c) . (p2 - p1))(p1 - c) + 
                    ((p1 - c) . (p1 - p2))(p2 - c)]];
    (p1 + p2)/2 + ((c - p1) . (-c + p2))/(2(-c + p2) . h)h
    ]

    (* Transformations between various models  *)
kleinToPoincare[x_List] := N[x/(1 + Sqrt[1 - x . x])] 
poincareToKlein[x_List] := N[(2x)/(1 + x . x)]
halfPlaneToPoincare[x_List] := 
  N[Composition[{Re[#], Im[#]}&, ((I - #)/(I + #))&, (# . {1, I})&][x]]
poincareToHalfPlane[x_List] := 
  N[Composition[{Re[#],Im[#]}&, I*((1 - #)/(1 + #))&, (# . {1, I})&][x]]

     (* LPoints in various models *)
PDPoint[{x_, y_}] :=
  LPoint[poincareToKlein[{x, y}]]/; (N[x^2 + y^2] < 1)
PDPoint[x___] := Message[PDPoint::badArg]
     (* Internal model for LPoint is Klein Disk *)
KDPoint[{x_, y_}] := 
  LPoint[{x, y}] /; (N[x^2 + y^2] < 1)              
KDPoint[x___] := Message[KDPoint::badArg]
HPPoint[{x_, y_}] := 
  LPoint[poincareToKlein[halfPlaneToPoincare[{x, y}]]] /; (N[y] > 0)
HPPoint[x___] := Message[HPPoint::badArg]

                   (* ISOMETRIC TRANSFORMATIONS *)
              (* (preformed in internal KleinDisk model *)

L2Reflection[a_LPoint][x_LPoint] :=
  Module[{a1, mat},
    a1 = LPointToH[a];
    mat = 
      IdentityMatrix[3] -
      (2/minDot[a1, a1])Map[List,a1] . {a1} . DiagonalMatrix[{1, 1, -1}];
    hToLPoint[mat . LPointToH[x]]
  ]

L2Reflection[LLine[{a_LPoint, b_LPoint}]][x_LPoint] :=
  Module[{ab, x1 = x[[1]], y1 = a[[1]] - b[[1]]},
    ab = Cross[LPointToH[a], LPointToH[b]];
    If[ab[[3]] == 0,
      (* If ab contains the origin the reflection is euclidean *)
      LPoint[((2 x1 . y1)/(y1 . y1))y1 - x1],
      L2Reflection[
         hToLPoint[Cross[ab, Cross[LPointToH[x], ab {1, 1, -1}]]]
      ][x]
    ]]

L2Translation[a_LPoint, b_LPoint][c_LPoint] :=
  Module[{ha, hb, mid},
    ha = LPointToH[a];
    hb = LPointToH[b];
    mid = 
      N[ha Sqrt[minDot[hb, hb] minDot[ha, hb]] + 
        hb Sqrt[minDot[ha, ha] minDot[ha, hb]]];
    L2Reflection[hToLPoint[mid]][L2Reflection[a][c]]]

   (* translate to the origin, do the euclidean rotation 
      and translate back *)
L2Rotation[a_LPoint, \[Theta]_][x_LPoint] :=                 
  Composition[
      L2Translation[LPoint[{0, 0}], a][#] &,
      LPoint[{{Cos[\[Theta]], Sin[\[Theta]]},{-Sin[\[Theta]],Cos[\[Theta]]}} . (#[[1]])] &,
      L2Translation[a, LPoint[{0, 0}]][#] &
      ][x]

l2Isometries = {L2Reflection, L2Translation, L2Rotation}

    (* LPrimitives are mapped in the same way by l2Isometries *) 
LCircle /: (f_)[y___][LCircle[x_LPoint, r_]] := 
  LCircle[f[y][x], r] /; MemberQ[l2Isometries, f]
LDisk /: (f_)[y___][LDisk[x_LPoint, r_]] := 
  LDisk[f[y][x], r] /; MemberQ[l2Isometries, f]
LLine /: (f_)[y___][LLine[l_List]] := 
  LLine[(f[y][#] &) /@  l] /; MemberQ[l2Isometries, f]
LPolygon /: (f_)[y___][LPolygon[l_List]] := 
  LPolygon[(f[y][#] &) /@  l] /; MemberQ[l2Isometries, f] 

   (* coordinates of LPoint in the model given in opts *)
coords[LPoint[x_], opts___] :=  
  Switch[Model /. {opts} /. Options[LToGraphics],
    KleinDisk, x,
    PoincareDisk, kleinToPoincare[x],
    HalfPlane, poincareToHalfPlane[kleinToPoincare[x]]
    ]

              (* IMPLEMENTATION OF LToGraphics  *) 

LToGraphics[x_LPoint, opts___] := Point[coords[x, opts]]

   (* makeArc-returns list of vertices of a polygonal line 
      approximating LLine ab in the model given in opts *)
   (* to be used ONLY when model is PoincareDisk and HalfPlane !!! *)
makeArc[a_LPoint, b_LPoint, opts___] :=                             
  Module[{c, a1, b1, angle, r, e1, e2, model, plotPoints},
    model = Model /. {opts} /. Options[LToGraphics];
    plotPoints = PlotPoints /. {opts} /. Options[LToGraphics];
    a1 = coords[a, opts];
    b1 = coords[b, opts];
    If[((model === PoincareDisk) && (Abs[Det[{a1, b1}]] < err)) || 
       ((model === HalfPlane) && (Abs[a1[[1]] - b1[[1]]] < err)),  
    (* arc is almost straight or straight-return line segment *)
       {a1, b1},
    (* otherwise- calculate the centre, radius and angle of the arc *)
       c = CCC[a1,b1,If[model === PoincareDisk,a1/(a1 . a1),a1{1, -1}]];
       r = norm[c - a1];
       angle = ArcCos[((a1 - c) . (b1 - c))/(r^2)];
       e1 = normalize[a1 - c];
       e2 = normalize[r^2(b1 - c)-((b1 - c) . (a1 - c))(a1 - c)];
       Table[c + r*( e1 Cos[\[Phi]]+e2 Sin[\[Phi]]), 
             {\[Phi], 0, angle, angle/plotPoints}]                      
      ] (* If *)
   ]

LToGraphics[LLine[points_List], opts___] :=
  Module[{vertices, model, l = Length[points]},
      model = Model /. {opts} /. Options[LToGraphics];
      If[model === KleinDisk,  
      (* all lines in KleinDisk are segments *)                         
        vertices = (points /. {LPoint[x_] :> x}),
                                    
      (* otherwise-draw arc by arc *)
        vertices = makeArc[points[[1]], points[[2]], opts];                 
        For[i = 2, i <= l - 1, i++,
          vertices =
             Join[vertices, 
               Drop[makeArc[points[[i]], points[[i + 1]], opts], 1]]    
             ] 
         ];
       Line[vertices]
   ] /; (Length[points] > 1)

LToGraphics[LPolygon[points_List], opts___] :=  
  Module[{vertices},  
      vertices = LToGraphics[LLine[points], opts][[1]];
      If[First[vertices] == Last[vertices],
     (* if polygon is closed - return as it is *) 
        Polygon[vertices],
     (* otherwise - close the polygon *)   
        Polygon[
           Join[vertices,
              Drop[LToGraphics[
                   LLine[{Last[points], First[points]}], opts][[1]],1]
           ] (* Join *)          
         ] (* Polygon *)
      ] (* If *)
  ] /; (Length[points] > 2)

LToGraphics[LCircle[lCentre_LPoint, r_], opts___] :=         
  Module[{model, plotPoints, dist, eCentre,rr, r1, r2, e1, e2},
       model = Model /. {opts} /.Options[LToGraphics];  
       If[model === PoincareDisk, rr = r/2, rr = r];
       c = coords[lCentre, Model -> model];
       dist = (norm[c]/(Cosh[rr]^2 - norm[c]^2  Sinh[rr]^2));
     (* Euclidean centre of a circle (ellipse) *)  
       eCentre = 
          If[model === HalfPlane, c {1, Cosh[rr]}, normalize[c] dist ];
     (* Euclidean (shorter) radius  of a circle (ellipse) *)    
       r1 = If[model === HalfPlane,
               c[[2]] Sinh[rr],
               1/2 dist(1-norm[c]^2)Sinh[2 rr]]; 
     (* LCircle is an ellipse in KleinDisk, a circle otherwise *)            
       Switch[model,
              HalfPlane, Circle[eCentre, r1],                             
              PoincareDisk, Circle[eCentre, r1],
              KleinDisk,        
                  e1 = normalize[c];
                  e2 = {-e1[[2]], e1[[1]]};
                  r2 = Sqrt[1 - dist^2] Tanh[rr];
                  plotPoints = 
                     PlotPoints /. {opts} /. Options[LToGraphics]; 
                  Line[Table[eCentre+r1 e1 Cos[\[Phi]]+r2 e2 Sin[\[Phi]],
                             {\[Phi], 0, 2 Pi, (2 Pi)/plotPoints}]]
        ] (* switch *)
   ] (* module *)

   (* LDisk is just a Polygon version of LCircle *)
LToGraphics[LDisk[lCentre_LPoint, r_], opts___] := 
  Module[{circle, plotPoints},
       circle = LToGraphics[LCircle[lCentre, r], opts];
       plotPoints = PlotPoints /. {opts} /. Options[LToGraphics];
       Switch[Model /. {opts} /. Options[LToGraphics],
              HalfPlane, Disk @@ circle,
              PoincareDisk, Disk @@ circle,
              KleinDisk, Polygon @@ circle
       ]
  ]  

 LToGraphics[l_List, opts___] := Map[LToGraphics[#, opts]&, l]        

End[]      
EndPackage[]
