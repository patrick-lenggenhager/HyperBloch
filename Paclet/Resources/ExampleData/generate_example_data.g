# load the HyperCells package (https://github.com/patrick-lenggenhager/HyperCells)
LoadPackage( "HyperCells" );

# {8,8} lattice

# set up (proper) triangle group
tg := ProperTriangleGroup( [ 2, 8, 8 ] );

# specify the quotient defining the primitive cell
qpc := TGQuotient( [ 2, 6 ] );

# construct symmetric primitive cell
cgpc := TGCellGraph( tg, qpc, 3 : simplify := 5 );
Export( cgpc, "(2,8,8)_T2.6_3.hcc" ); # export

# elementary nearest-neighbor model
model := TessellationModelGraph( cgpc, true : simplify := 5 );
Export( model, "{8,8}-tess_T2.6_3.hcm" ); # export

# construct symmetric supercell
sc := TGCellSymmetric( tg, TGQuotient( [ 3, 11 ] ), 3 );

# extend the model defined on the primitive cell to the supercell
scmodel := TGSuperCellModelGraph( model, sc : simplify := 0 );
Export( scmodel, "{8,8}-tess_T2.6_3_sc-T3.11.hcs" ); # export

# {8,3} lattice

# set up (proper) triangle group
tg := ProperTriangleGroup( [ 2, 3, 8 ] );

# specify the quotient defining the primitive cell
qpc := TGQuotient( [ 2, 1 ] );

# construct symmetric primitive cell
cgpc := TGCellGraph( tg, qpc, 3 : simplify := 5 );
Export( cgpc, "(2,3,8)_T2.1_3.hcc" ); # export

# elementary nearest-neighbor model
model := TessellationModelGraph( cgpc, false : simplify := 5 );
Export( model, "{8,3}-tess_T2.1_3.hcm" ); # export

# construct symmetric supercell
sc := TGCellSymmetric( tg, TGQuotient( [ 5, 1 ] ), 3 );

# extend the model defined on the primitive cell to the supercell
scmodel := TGSuperCellModelGraph( model, sc : simplify := 0 );
Export( scmodel, "{8,3}-tess_T2.1_3_sc-T5.1.hcs" ); # export