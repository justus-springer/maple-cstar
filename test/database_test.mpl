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
with(LinearAlgebra):

# Hard coded indices for the CSV table

$define DB_ID 1
$define DB_FORMAT 2
$define DB_P 3
$define DB_Q 5
$define DB_ANTICAN 8
$define DB_PICARDNUMBER 14

local CSVData, t, row, r, P, expected_Q, evaluated_Q, i, j, expected_antican, d, Sigma:

CSVData := ImportMatrix("gorenstein_database.txt", delimiter = ";"):

t := time():
# We skip the first row, as these are just labels
for i from 2 to RowDimension(CSVData) do
    row := CSVData[i]:
    r := nops(parse(row[DB_FORMAT])) - 1:
    P := PMatrix(r, Matrix(parse(row[DB_P]))):

    # Test if the degree matrix is as expected
    expected_Q := parse(row[DB_Q]):
    evaluated_Q := map(x -> convert(x, list), [Row(getQ(P), [seq(1 .. RowDimension(getQ(P)))])]):
    Test(evaluated_Q, expected_Q, label = cat("Variety id: ", row[DB_ID], ", Degree matrix test")):

    # Test if the anticanonical class is as expected
    expected_antican := parse(row[DB_ANTICAN]):
    Test(convert(getAnticanClass(P), list), expected_antican, label = cat("Variety id: ", row[DB_ID], ", Antican test")):
    
    # Gorenstein test
    d := P:-n + P:-m;
    Sigma := convert(combinat[powerset]({seq(1..d)}) minus {{seq(1..d)}}, list):
    X := TVarOne(P, Sigma);
    Test(isGorenstein(X), true, label = cat("Variety id: ", row[DB_ID], ", Gorenstein test")):

end do:

t := time() - t:
print(cat("Elapsed time = ", t, "s")):