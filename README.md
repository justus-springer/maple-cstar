# complexity-one

A **[Maple](https://www.maplesoft.com/products/maple/)** package for rational
T-Varieties with a torus action of complexity one. A special emphasis is on the
two-dimensional case (C*-surfaces). It makes use of the **[convex package](https://www.math.uwo.ca/faculty/franz/convex/)** by Matthias Franz.

## Mathematical background

This software package is based on a combinatorial approach to complexity one
varieties that is described in Huggenberger's thesis **[Fano Varieties with a
Torus Action of Complexity
One](https://publikationen.uni-tuebingen.de/xmlui/handle/10900/49921)**.
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
  self-intersection.

## Installation

Before installing this package, make sure the *[convex package](https://www.math.uwo.ca/faculty/franz/convex/)* is properly installed by typing `with(convex);` into a maple prompt.

To install the package, download the Maple Library Archive file from the
**[Releases section](https://github.com/justus-springer/complexity-one/releases)**
and put it into a location present in your Maple's `libname` variable. To see
the current value of `libname`, just type `libname;` into a maple prompt and
hit enter. For more information, see the associated Maple **[help page on
libname](https://www.maplesoft.com/support/help/Maple/view.aspx?path=libname)**.
In particular, note that you cak automatically asign a custom value to
`libname` at startup by editing **[Maple's initialization
file](https://www.maplesoft.com/support/help/Maple/view.aspx?path=worksheet%2freference%2finitialization)** (For instance, you can add a custom path where you can put all your user-added libraries).

Alternatively, you can also download the source code of this repository. Starting a Maple session in the directory where the code is, you can run

```
read "ComplexityOne.mpl";
```

to read in the source code of the package and start using it immediately.

## Usage

Type `with(ComplexityOnePackage);` into a maple prompt to make the procedures of the package available. 

