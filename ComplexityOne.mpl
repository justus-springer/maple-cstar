ComplexityOne := module()

option package;

export 
    PFormat, 
    AdmissibleOperation, 
    PMatrix, 
    ComplexityOneVariety, 
    
    # Database stuff
    FindInDatabase, 
    ImportComplexityOneVarietyList, 
    ImportComplexityOneVariety, 
    ExportComplexityOneVarietyList, 
    performOnDatabase;

local
    applyPermToList,
    invariantPermutations,
    sortLexComparison,
    sortLex;

uses LinearAlgebra, Database[SQLite];

$include "Tools.mm"
$include "PFormat.mm"
$include "AdmissibleOperation.mm"
$include "PMatrix.mm"
$include "ComplexityOneVariety.mm"
$include "DatabaseTools.mm"

end module:
