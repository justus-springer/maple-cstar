# complexity-one

A **[Maple](https://www.maplesoft.com/products/maple/)** package for rational
T-Varieties with a torus action of complexity one. A special emphasis is on the
two-dimensional case (C*-surfaces). It makes use of the **[convex package](https://www.math.uwo.ca/faculty/franz/convex/)** by Matthias Franz.

## Mathematical background

This software package is based on a combinatorial approach to complexity one
varieties that is described in **[Huggenberger's thesis]
(https://publikationen.uni-tuebingen.de/xmlui/handle/10900/49921)**.
Another reference is Chapter 1.5 of **[Nicolussi's
thesis](https://publikationen.uni-tuebingen.de/xmlui/handle/10900/77235)**. For
background on C*-surfaces and the isomorphy problem, consult for instance
**[arxiv:2302.03095](https://arxiv.org/abs/2302.03095)**.

## Features

- Computation of divisor class group, canonical divisor class, effective cone,
  moving cone and ample cone of divisor classes,
- Tests for being Fano, Q-Factorial and Q-Gorenstein, computation of Gorenstein
  Index and Picard Index,
- Tests for whether or not two sets of defining data describe isomorphic
  complexity-one varieties,
- For C*-surfaces: Efficient computation of a normal form for the defining data,
- For C*-surfaces: Computation of intersection numbers and anticanonical
  self-intersection, test for log terminality, determining the singularity type, computing a resolution of singularities.

## Installation

Before installing this package, make sure the **[convex package](https://www.math.uwo.ca/faculty/franz/convex/)** is properly installed by typing `with(convex);` into a maple prompt.

To install the package, download the Maple Library Archive file from the
**[Releases section](https://github.com/justus-springer/complexity-one/releases)**
and put it into a location present in your Maple's `libname` variable. To see
the current value of `libname`, just type `libname;` into a maple prompt and
hit enter. For more information, see the associated Maple **[help page on
libname](https://www.maplesoft.com/support/help/Maple/view.aspx?path=libname)**.
In particular, note that you can automatically asign a custom value to
`libname` at startup by editing **[Maple's initialization
file](https://www.maplesoft.com/support/help/Maple/view.aspx?path=worksheet%2freference%2finitialization)** (For instance, you can add a custom path where you can put all your user-added libraries).

Alternatively, you can use the following commands to download the source code from this repository and read it directly into maple: 

```
git clone https://github.com/justus-springer/complexity-one
cd complexity-one
maple
```

Now, from within maple, run

```
read "ComplexityOne.mpl";
```

## Quick start

Below is a very quick start guide. For more information, see the file **[tutorial.txt](https://github.com/justus-springer/complexity-one/blob/main/tutorial.txt)**
in this repository as well as the documentation in the source code.

First, load `convex` and `ComplexityOnePackage`.

```
with(convex);
with(ComplexityOnePackage);
```

Define a Complexity-one variety from a P-Matrix (in this case, a C*-surface)

```
P := PMatrix(1, <-1, -5, 2, 0; -1, -5, 0, 2; 0, -6, 1, 1>);
X := ComplexityOneVariety(P);
```

Some invariants are immediately accessible from `P` and `X`:

```
P:-dim;
P:-classGroupRank;
P:-case;
X:-Sigma;
```

Others can be computed (and are then cached) from various `get` functions:

```
getClassGroup(P);
getCanonicalDivisorClass(P);
isLogTerminal(P);
getSingularityType(P);
isQfactorial(X);
isQgorenstein(X);
isGorenstein(X);
isFano(X);
getGorensteinIndex(X);
getPicardIndex(X);
getIntersectionMatrix(X);
getAnticanonicalSelfIntersection(X);
```

Getting a minimal resolution of singularities and checking that it is factorial:

```
Xres := minimalResolution(X);
isFactorial(Xres);
Xres:-P:-mat;
```
