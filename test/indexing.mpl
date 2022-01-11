$include "../complexity_one.mpl"

with(ComplexityOne):
with(CodeTools):

(*
Some tests about `singleToDoubleIndex` and `doubleToSingleIndex`
*)

f := PFormat([2,3,1,2],3,2);

Test([singleToDoubleIndex(f,1)], [1,1], quiet, label = "singleToDoubleIndex-1");
Test([singleToDoubleIndex(f,2)], [1,2], quiet, label = "singleToDoubleIndex-2");
Test([singleToDoubleIndex(f,3)], [2,1], quiet, label = "singleToDoubleIndex-3");
Test([singleToDoubleIndex(f,4)], [2,2], quiet, label = "singleToDoubleIndex-4");
Test([singleToDoubleIndex(f,5)], [2,3], quiet, label = "singleToDoubleIndex-5");
Test([singleToDoubleIndex(f,6)], [3,1], quiet, label = "singleToDoubleIndex-6");
Test([singleToDoubleIndex(f,7)], [4,1], quiet, label = "singleToDoubleIndex-7");
Test([singleToDoubleIndex(f,8)], [4,2], quiet, label = "singleToDoubleIndex-8");
Test([singleToDoubleIndex(f,9)], [-1,1], quiet, label = "singleToDoubleIndex-9");
Test([singleToDoubleIndex(f,10)], [-1,2], quiet, label = "singleToDoubleIndex-10");
Test([singleToDoubleIndex(f,11)], [-1,3], quiet, label = "singleToDoubleIndex-11");

Test([singleToDoubleIndex(f,0)], "", 'testerror', quiet, label = "singleToDoubleIndex-12");
Test([singleToDoubleIndex(f,-3)], "", 'testerror', quiet, label = "singleToDoubleIndex-13");
Test([singleToDoubleIndex(f,12)], "", 'testerror', quiet, label = "singleToDoubleIndex-14");

Test(doubleToSingleIndex(f,1,1), 1, quiet, label = "doubleToSingleIndex-1");
Test(doubleToSingleIndex(f,1,2), 2, quiet, label = "doubleToSingleIndex-2");
Test(doubleToSingleIndex(f,2,1), 3, quiet, label = "doubleToSingleIndex-3");
Test(doubleToSingleIndex(f,2,2), 4, quiet, label = "doubleToSingleIndex-4");
Test(doubleToSingleIndex(f,2,3), 5, quiet, label = "doubleToSingleIndex-5");
Test(doubleToSingleIndex(f,3,1), 6, quiet, label = "doubleToSingleIndex-6");
Test(doubleToSingleIndex(f,4,1), 7, quiet, label = "doubleToSingleIndex-7");
Test(doubleToSingleIndex(f,4,2), 8, quiet, label = "doubleToSingleIndex-8");
Test(doubleToSingleIndex(f,-1,1), 9, quiet, label = "doubleToSingleIndex-9");
Test(doubleToSingleIndex(f,-1,2), 10, quiet, label = "doubleToSingleIndex-10");
Test(doubleToSingleIndex(f,-1,3), 11, quiet, label = "doubleToSingleIndex-11");

Test(doubleToSingleIndex(f,1,3), "", 'testerror', quiet, label = "doubleToSingleIndex-12");
Test(doubleToSingleIndex(f,3,2), "", 'testerror', quiet, label = "doubleToSingleIndex-13");
Test(doubleToSingleIndex(f,1,0), "", 'testerror', quiet, label = "doubleToSingleIndex-14");
Test(doubleToSingleIndex(f,-1,0), "", 'testerror', quiet, label = "doubleToSingleIndex-15");
Test(doubleToSingleIndex(f,-1,4), "", 'testerror', quiet, label = "doubleToSingleIndex-16");