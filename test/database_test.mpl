(*
This file reads in all of the P-matrices in `gorenstein_database.csv` and verifies that all
of the entries are correct.

Currently, we only test if the degree matrix, anticanonical class and the Picard number 
is as expected, and if they all are gorenstein.
*)

$include "../complexity_one.mpl"

with(convex):
with(ComplexityOne):
with(CodeTools):

local t, Ps, i, P, numColumns, Sigma, X:

t := time():

printf("Parsing and verifying gorenstein_database.txt...\n"):

Ps := parseCSV("../gorenstein_database.txt"):

printf("Parsing succesfull. Performing gorenstein tests...\n"):

# Check the gorenstein condition for every P-matrix.
for i from 1 to nops(Ps) do

    P := Ps[i];
    numColumns := P:-n + P:-m;
    Sigma := convert(combinat[powerset]({seq(1..numColumns)}) minus {{seq(1..numColumns)}}, list):
    X := TVarOne(P, Sigma);
    Test(isGorenstein(X), true, label = cat("Row: ", i, ", Gorenstein test")):

end do:

t := time() - t:
printf("Elapsed time = %gs\n", t):