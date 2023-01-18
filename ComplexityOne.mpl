ComplexityOnePackage := module()

option package;

export 
    PFormat, 
    AdmissibleOperation, 
    PMatrix, 
    ComplexityOneVariety, 

    # Tools
    swapCase,
    platonicityType,
    canonicalBasisVector,
    applyPermToList,
    invariantPermutations,
    sortLexComparison,
    sortLex,
    intersectLines2D,
    hilbertBasisCone2D,
    getSingleDiscrepancy,
    imageFactorGroup,
    indexOfImage,
    dualHermiteForm,
    integerKernel,
    integerIntersectionBasis,
    integerIntersectionBasisList,
    
    # Database stuff
    CreateDatabase,
    ImportFromDatabase,
    ImportFromDatabaseID,
    FindInDatabase, 
    ExportToDatabase,
    UpdateInDatabase,
    performOnDatabase;

uses LinearAlgebra, Database[SQLite];

$include "Tools.mm"
$include "PFormat.mm"
$include "AdmissibleOperation.mm"
$include "PMatrix.mm"
$include "ComplexityOneVariety.mm"
$include "DatabaseTools.mm"

end module:
