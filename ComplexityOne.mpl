ComplexityOnePackage := module()

option package;

export 
    PFormat, 
    AdmissibleOperation, 
    PMatrix, 
    ComplexityOneVariety, 

    # Tools
    swapCase,
    canonicalBasisVector,
    applyPermToList,
    invariantPermutations,
    sortLexComparison,
    sortLex,
    imageFactorGroup,
    indexOfImage,
    
    # Database stuff
    FindInDatabase, 
    ImportComplexityOneVarietyList, 
    ImportComplexityOneVariety, 
    ExportComplexityOneVarietyList, 
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
