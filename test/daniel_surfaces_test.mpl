(*
This file reads in all of the P-matrices in `databases/daniel_surfaces.txt` and verifies that
exactly the specified ones are gorenstein.
*)

$include "../complexity_one.mpl"

with(convex):
with(ComplexityOne):
with(CodeTools):

t1 := time():

printf("Parsing and verifying gorenstein_database.txt...\n"):

Ps := ImportPMatrixList("../databases/daniel_surfaces.txt"):
Xs := map(P -> TVarOne(P), Ps):

t2 := time():

# Exactly these indices correspond to gorenstein matrices
gorenstein_indices := 
[1, 10, 31, 36, 38, 39, 40, 44, 45, 51, 69, 83, 91, 92, 96,
 100, 101, 110, 126, 149, 150, 310, 311, 314, 326, 327, 332,
 347, 410, 437, 1017, 1020, 1033, 1055]:


printf("Parsing succesfull. Performing gorenstein tests...\n"):

for i from 1 to nops(Ps) do
    if i in gorenstein_indices then
        res := true;
    else
        res := false;
    end if;
    Test(isGorenstein(TVarOne(Ps[i])), res, quiet, label = cat("Row: ", i, ", Gorenstein test")):
end do:

t3 := time():
printf("Time for Parsing: %gs\n", t2 - t1):
printf("Time for gorentein tests: %gs\n", t3 - t2):