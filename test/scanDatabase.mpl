$include "../complexity_one.mpl"

with(ComplexityOne):
with(CodeTools):
with(LinearAlgebra):
with(convex):
with(Database[SQLite]):

db := Open("../surfaces.db");
tableName := "surfaces";

testVarietyConsistency := proc(X :: TVarOne, i :: integer)
    print(cat("Testing variety no. ", i));
    # Test class Group
    cachedClassGroup := getClassGroup(X:-P);
    computedClassGroup := getClassGroup(X:-P, 'forceCompute');
    Test(computedClassGroup, cachedClassGroup, quiet, label = cat("TEST_classGroup_", i));
    # Test degree matrix
    cachedDegreeMatrix := getQ(X:-P);
    computedDegreeMatrix := getQ(X:-P, 'forceCompute');
    Test(Equal(computedDegreeMatrix, cachedDegreeMatrix), true, quiet, label = cat("TEST_degreeMatrix_", i));
    # Test anticanonical class
    cachedAnticanClass := getAnticanClass(X:-P);
    computedAnticanClass := getAnticanClass(X:-P, 'forceCompute');
    Test(Equal(computedAnticanClass, cachedAnticanClass), true, quiet, label = cat("TEST_anticanClass_", i));
    # Test gorenstein index
    cachedGorensteinIndex := getGorensteinIndex(X);
    computedGorensteinIndex := getGorensteinIndex(X, 'forceCompute');
    Test(computedGorensteinIndex, cachedGorensteinIndex, quiet, label = cat("TEST_gorensteinIndex_", i));
    # Check if the variety occurs only once in the database
    Test(FindInDatabase(db, tableName, X), [i], quiet, quiet, label = cat("TEST_UNIQUE_", i));
end proc;

performOnDatabase(db, tableName, testVarietyConsistency);