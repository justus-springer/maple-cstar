(*
Some tests for `isBigCone` and `isLeafCone`
*)

with(ComplexityOnePackage):
with(CodeTools):

f1 := PFormat([2,1,1], 1, 1):
f2 := PFormat([3,2,2,1], 2, 1):

Test(isBigCone(f1, {1,3,4}), true, label = "isBigCone-1");
Test(isBigCone(f1, {2,3,4}), true, label = "isBigCone-2");
Test(isBigCone(f1, {1,3,4,5}), true, label = "isBigCone-3");
Test(isBigCone(f2, {1,4,6,8}), true, label = "isBigCone-4");
Test(isBigCone(f2, {2,5,7,8}), true, label = "isBigCone-5");
Test(isBigCone(f2, {2,3,4,6,7,8,10}), true, label = "isBigCone-6");
Test(isBigCone(f2, {1,5,6,8,9}), true, label = "isBigCone-7");

Test(isBigCone(f1, {3,4}), false, label = "isBigCone-8");
Test(isBigCone(f1, {1,2}), false,  quiet, label = "isBigCone-9");
Test(isBigCone(f1, {1,3,5}), false, label = "isBigCone-10");
Test(isBigCone(f1, {3,4,5}), false, label = "isBigCone-11");
Test(isBigCone(f2, {4,5,6,7,8,9,10}), false, label = "isBigCone-12");
Test(isBigCone(f2, {1,3,6,8}), false, label = "isBigCone-13");
Test(isBigCone(f2, {2,5,6,9}), false, label = "isBigCone-14");

Test(isLeafCone(f1, {1,2,5}), true, label = "isLeafCone-1");
Test(isLeafCone(f1, {3,5}), true, label = "isLeafCone-2");
Test(isLeafCone(f1, {3,4}), false, label = "isLeafCone-3");
Test(isLeafCone(f1, {1,3,5}), false, label = "isLeafCone-4");

Test(isLeafCone(f2, {1,3,9,10}), true, label = "isLeafCone-5");
Test(isLeafCone(f2, {8,9}), true, label = "isLeafCone-6");
Test(isLeafCone(f2, {1,4,9,10}), false, label = "isLeafCone-7");
Test(isLeafCone(f2, {7,8,9,10}), false, label = "isLeafCone-8");
