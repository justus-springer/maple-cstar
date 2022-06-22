with(ComplexityOnePackage):
with(CodeTools):

(*
Some tests about `singleToDoubleIndex` and `doubleToSingleIndex`
*)

f := PFormat([2,3,1,2],3,2);

Test([singleToDoubleIndex(f,1)], [0,1], label = "singleToDoubleIndex-1");
Test([singleToDoubleIndex(f,2)], [0,2], label = "singleToDoubleIndex-2");
Test([singleToDoubleIndex(f,3)], [1,1], label = "singleToDoubleIndex-3");
Test([singleToDoubleIndex(f,4)], [1,2], label = "singleToDoubleIndex-4");
Test([singleToDoubleIndex(f,5)], [1,3], label = "singleToDoubleIndex-5");
Test([singleToDoubleIndex(f,6)], [2,1], label = "singleToDoubleIndex-6");
Test([singleToDoubleIndex(f,7)], [3,1], label = "singleToDoubleIndex-7");
Test([singleToDoubleIndex(f,8)], [3,2], label = "singleToDoubleIndex-8");
Test([singleToDoubleIndex(f,9)], [-1,1], label = "singleToDoubleIndex-9");
Test([singleToDoubleIndex(f,10)], [-1,2], label = "singleToDoubleIndex-10");
Test([singleToDoubleIndex(f,11)], [-1,3], label = "singleToDoubleIndex-11");

Test([singleToDoubleIndex(f,0)], "", 'testerror', label = "singleToDoubleIndex-12");
Test([singleToDoubleIndex(f,-3)], "", 'testerror', label = "singleToDoubleIndex-13");
Test([singleToDoubleIndex(f,12)], "", 'testerror', label = "singleToDoubleIndex-14");

Test(doubleToSingleIndex(f,0,1), 1, label = "doubleToSingleIndex-1");
Test(doubleToSingleIndex(f,0,2), 2, label = "doubleToSingleIndex-2");
Test(doubleToSingleIndex(f,1,1), 3, label = "doubleToSingleIndex-3");
Test(doubleToSingleIndex(f,1,2), 4, label = "doubleToSingleIndex-4");
Test(doubleToSingleIndex(f,1,3), 5, label = "doubleToSingleIndex-5");
Test(doubleToSingleIndex(f,2,1), 6, label = "doubleToSingleIndex-6");
Test(doubleToSingleIndex(f,3,1), 7, label = "doubleToSingleIndex-7");
Test(doubleToSingleIndex(f,3,2), 8, label = "doubleToSingleIndex-8");
Test(doubleToSingleIndex(f,-1,1), 9, label = "doubleToSingleIndex-9");
Test(doubleToSingleIndex(f,-1,2), 10, label = "doubleToSingleIndex-10");
Test(doubleToSingleIndex(f,-1,3), 11, label = "doubleToSingleIndex-11");

Test(doubleToSingleIndex(f,0,3), "", 'testerror', label = "doubleToSingleIndex-12");
Test(doubleToSingleIndex(f,2,2), "", 'testerror', label = "doubleToSingleIndex-13");
Test(doubleToSingleIndex(f,0,0), "", 'testerror', label = "doubleToSingleIndex-14");
Test(doubleToSingleIndex(f,-1,0), "", 'testerror', label = "doubleToSingleIndex-15");
Test(doubleToSingleIndex(f,-1,4), "", 'testerror', label = "doubleToSingleIndex-16");