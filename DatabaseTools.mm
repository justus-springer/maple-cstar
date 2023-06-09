(*
Some tools for handling a database of K*-surfaces.
*)

(*
Creates a SQLite database for storing K*-surfaces with a given filepath
*)
CreateDatabase := proc(filePath :: string, tableName :: string)
    local db;
    db := Open(filePath);
    Execute(db, cat("CREATE TABLE ", tableName, "(r, ns, n, m, s, lss, P, dimension, classGroupRank, classGroup, degreeMatrix, canonicalDivisorClass, ambientFan, maximalXCones, gorensteinIndex, isGorenstein, orderedLss, effectiveCone, movingCone, ampleCone, isFano, variables, monomials, relations, intersectionMatrix, anticanonicalSelfIntersectionFraction, anticanonicalSelfIntersectionFloat, isToric, isIrredundant, isQfactorial, isFactorial, localPicardIndices, picardIndex, isQgorenstein, localGorensteinIndices, localGorensteinQuotients, gorensteinQuotient, case_, orientation, isLogTerminal, singularityType);"));
    return db;
end proc;

(*
Imports a list of K*-surfaces from a given SQLite statement
*)
ImportFromDatabase := proc(stmt)
    local columns, INDEX_S, INDEX_P, i, clm, result;

    # Find the columns indices for `P` and `s`
    columns := ColumnNames(stmt);
    for i from 0 to nops(columns)-1 do
        clm := columns[i+1];
        if clm = "s" then
            INDEX_S := i;
        elif clm = "P" then
            INDEX_P := i;
        end if;
    end do;

    if not type(INDEX_S, integer) or not type(INDEX_P, integer) then
        error "Not enough data. You must provide the columns `P` and `s`.";
    end if;

    result := [];
    while Step(stmt) = RESULT_ROW do
        result := [op(result), ComplexityOneVariety(PMatrix(Fetch(stmt, INDEX_S), Matrix(parse(Fetch(stmt, INDEX_P))), 'skipChecks'))];
    end do;

    return result;

end proc;

(*
Imports a single K*-surface from a database with a given rowid.
To import multiple varieties using an arbitrary SQLite query, use `ImportComplexityOneVarietyList`
*)
ImportFromDatabaseID := proc(db, tableName :: string, rowid :: integer)
    ImportFromDatabase(Prepare(db, cat("SELECT * FROM ", tableName, " WHERE rowid = ", rowid)))[1];
end proc;

FindInDatabase := proc(db, tableName :: string, X0 :: ComplexityOneVariety)
    local X, stmt;

    if not X0:-P:-dim = 2 then
        error "This method only works for K*-surfaces."
    end if;

    X := normalForm(removeErasableBlocks(X0));

    stmt := Prepare(db, 
        cat("SELECT rowid FROM ", tableName, " WHERE ",
        "P = \"", convert(X:-P:-mat, list, nested), "\""));

    return convert(FetchAll(stmt), list);

end proc;

(*
Given a list of varieties of complexity one `Xs` and a SQLite database connection `db`, this function
inserts those K*-surfaces from `Xs` into the database that are not already present.
*)
ExportToDatabase := proc(connection, tableName :: string, Xs :: list(ComplexityOneVariety))
    local k, X, P, stmt, i, knownCount, numberOfColumns;
    knownCount := 0;

    for k from 1 to nops(Xs) do
        X := normalForm(removeErasableBlocks(Xs[k]));
        
        # Only insert `X` if it is not already present.
        # Can be disabled by the 'noChecks' option
        if 'noChecks' in [_passed] or FindInDatabase(connection, tableName, X) = [] then
            P := X:-P;
            numberOfColumns := if isToric(P) then 41 else 42 end if;
            stmt := Prepare(connection, cat("INSERT INTO ", tableName ," VALUES (", cat("?," $ numberOfColumns - 1), "?)"));
            Bind(stmt, 1, P:-r, valuetype = "integer");
            Bind(stmt, 2, convert(convert(P:-ns, list), string), valuetype = "text");
            Bind(stmt, 3, P:-n, valuetype = "integer");
            Bind(stmt, 4, P:-m, valuetype = "integer");
            Bind(stmt, 5, P:-s, valuetype = "integer");
            Bind(stmt, 6, convert(convert(P:-lss, list, nested), string), valuetype = "text");
            Bind(stmt, 7, convert([seq(convert(Row(P:-mat, i), list), i = 1 .. RowDimension(P:-mat))], string), valuetype = "text");
            Bind(stmt, 8, P:-dim, valuetype = "integer");
            Bind(stmt, 9, P:-classGroupRank, valuetype = "integer");
            Bind(stmt, 10, convert(getClassGroup(P), string), valuetype = "text");
            Bind(stmt, 11, convert([seq(convert(Row(getDegreeMatrix(P), i), list), i = 1 .. RowDimension(getDegreeMatrix(P)))], string), valuetype = "text");
            Bind(stmt, 12, convert(convert(getCanonicalDivisorClass(P), list), string), valuetype = "text");
            Bind(stmt, 13, convert(X:-Sigma, string), valuetype = "text");
            Bind(stmt, 14, convert(getMaximalXCones(X), string), valuetype = "text");
            Bind(stmt, 15, getGorensteinIndex(X), valuetype = "integer");
            Bind(stmt, 16, isGorenstein(X), valuetype = "integer");
            Bind(stmt, 17, convert(convert(sortColumnsByLss(P):-lss, list, nested), string), valuetype = "text");
            Bind(stmt, 18, convert(rays(getEffectiveCone(P)), string), valuetype = "text");
            Bind(stmt, 19, convert(rays(getMovingCone(P)), string), valuetype = "text");
            Bind(stmt, 20, convert(rays(getAmpleCone(X)), string), valuetype = "text");
            Bind(stmt, 21, isFano(X), valuetype = "integer");
            Bind(stmt, 22, convert(X:-variables, string), valuetype = "text");
            Bind(stmt, 23, convert(X:-monomials, string), valuetype = "text");
            Bind(stmt, 24, convert(X:-relations, string), valuetype = "text");
            Bind(stmt, 25, convert(convert(map(x -> [numer(x), denom(x)], getintersectionMatrix(X)), list, nested), string), valuetype = "text");
            Bind(stmt, 26, convert([numer(getAnticanonicalSelfIntersection(X)), denom(getAnticanonicalSelfIntersection(X))], string), valuetype = "text");
            Bind(stmt, 27, convert(getAnticanonicalSelfIntersection(X), hfloat), valuetype = "float");
            Bind(stmt, 28, isToric(P), valuetype = "integer");
            Bind(stmt, 29, isIrredundant(P), valuetype = "integer");
            Bind(stmt, 30, isQfactorial(X), valuetype = "integer");
            Bind(stmt, 31, isFactorial(X), valuetype = "integer");
            Bind(stmt, 32, convert(getLocalPicardIndices(X), string), valuetype = "text");
            Bind(stmt, 33, getPicardIndex(X), valuetype = "integer");
            Bind(stmt, 34, isQgorenstein(X), valuetype = "integer");
            Bind(stmt, 35, convert(getLocalGorensteinIndices(X), string), valuetype = "text");
            Bind(stmt, 36, convert(getLocalGorensteinQuotients(X), string), valuetype = "text");
            Bind(stmt, 37, getGorensteinQuotient(X), valuetype = "integer");
            Bind(stmt, 38, X:-P:-case, valuetype = "text");
            Bind(stmt, 39, X:-P:-orientation, valuetype = "integer");
            Bind(stmt, 40, isLogTerminal(X:-P), valuetype = "integer");
            if not isToric(P) then
                Bind(stmt, 41, getSingularityType(X:-P), valuetype = "text");
            end if;

            Step(stmt);
            Finalize(stmt);
            if 'logging' in [_passed] then
                print(cat("Variety no. ", k, " inserted."));
            end if;
        else
            knownCount := knownCount + 1;
            if 'logging' in [_passed] then
                print(cat("Variety no. ", k, " already known."));
            end if;
        end if;
    end do;

    if 'logging' in [_passed] then
        if 'noChecks' in [_passed] then
            print(cat("Written ", nops(Xs), "varieties to the database"));
        else
            print(cat("Supplied ", nops(Xs), " varieties. Already known: ", knownCount, ". Inserted ", nops(Xs) - knownCount, " new varieties into the database."));
        end if;
    end if;
    
end proc;

UpdateInDatabase := proc(connection, tableName :: string, rowid :: integer, X :: ComplexityOneVariety)
    local P, stmt;
    P := X:-P;    
    stmt := cat("UPDATE ", tableName, " SET ",
        "r = ", P:-r, ", ",
        "ns = \"", P:-ns, "\", ",
        "n = ", P:-n, ", ",
        "m = ", P:-m, ", ",
        "s = ", P:-s, ", ",
        "lss = \"", P:-lss, "\", ",
        "P = \"", convert(P:-mat, list, nested), "\", ",
        "dimension = ", P:-dim, ", ",
        "classGroupRank = ", P:-classGroupRank, ", ",
        "classGroup = \"", getClassGroup(P), "\", ",
        "degreeMatrix = \"", convert(getDegreeMatrix(P), list, nested), "\", ",
        "anticanClass = \"", convert(getCanonicalDivisorClass(P), list), "\", ",
        "ambientFan = \"", X:-Sigma, "\", ",
        "maximalXCones = \"", getMaximalXCones(X), "\", ",
        "gorensteinIndex = ", getGorensteinIndex(X), ", ",
        "isGorenstein = ", if isGorenstein(X) then 1 else 0 end if, ", ",
        "orderedLss = \"", sortColumnsByLss(P):-lss, "\", ",
        "effectiveCone = \"", rays(getEffectiveCone(P)), "\", ",
        "movingCone = \"", rays(getMovingCone(P)), "\", ",
        "ampleCone = \"", rays(getAmpleCone(X)), "\", ",
        "isFano = ", if isFano(X) then 1 else 0 end if, ", ",
        "variables = \"", X:-variables, "\", ",
        "monomials = \"", X:-monomials, "\", ",
        "relations = \"", X:-relations, "\", ",
        "intersectionMatrix = \"", convert(map(x -> [numer(x), denom(x)], getintersectionMatrix(X)), list, nested), "\", ",
        "anticanonicalSelfIntersectionFraction = \"", [numer(getAnticanonicalSelfIntersection(X)), denom(getAnticanonicalSelfIntersection(X))], "\", ",
        "anticanonicalSelfIntersectionFloat = ", convert(getAnticanonicalSelfIntersection(X), hfloat), ", ",
        "isToric = ", if isToric(P) then 1 else 0 end if, ", ",
        "isIrredundant = ", if isIrredundant(P) then 1 else 0 end if, ", ",
        "isQfactorial = ", if isQfactorial(X) then 1 else 0 end if, ", ",
        "isFactorial = ", if isFactorial(X) then 1 else 0 end if, ", ",
        "localPicardIndices = \"", getLocalPicardIndices(X), "\", ",
        "picardIndex = ", getPicardIndex(X), ", ",
        "isQgorenstein = ", if isQgorenstein(X) then 1 else 0 end if, ", ",
        "localGorensteinIndices = \"", getLocalGorensteinIndices(X), "\", ",
        "localGorensteinQuotients = \"", getLocalGorensteinQuotients(X), "\", ",
        "gorensteinQuotient = ", getGorensteinQuotient(X), ", ",
        "case_ = \"", P:-case, "\", ",
        "orientation = ", P:-orientation, " ",
        "isLogTerminal = ", if isLogTerminal(P) then 1 else 0 end if,
        "singularityType = \"", getSingularityType(P), "\" "
        "WHERE rowid = ", rowid
        );
    Execute(db, stmt);
end proc;

(*
Helper function to perform an operation on every variety in the database, without loading them into memory all at once
(which can be unpractical, when the database is large).
The parameter `f` is a function taking a `ComplexityOneVariety` as first argument and an integer as second argument, which is the
position the variety occurs in the database.
*)
performOnDatabase := proc(db, tableName :: string, f, step := 1000)
    local numberOfEntries, numberOfSteps, i, offset, Xs, j;
    # First, count the number of database entries
    numberOfEntries := FetchAll(Prepare(db, cat("SELECT COUNT(*) FROM ", tableName)))[1,1];
    numberOfSteps := ceil(numberOfEntries / step);
    for i from 1 to numberOfSteps do 
        offset := (i - 1) * step;
        Xs := ImportFromDatabase(Prepare(db, cat("SELECT * FROM ", tableName, " LIMIT ", step, " OFFSET ", offset)));
        for j from 1 to nops(Xs) do
            f(Xs[j], offset + j);
        end do;
    end do;
end proc;
