
# /progsys/maple/V11/bin/maple -q < standarddiffop.mpl

dirs := [[0, n],
         [ e, n ],
         [ e, 0 ],
         [ e, -s ],
         [ 0, -s ],
         [ -w, -s ],
         [ -w, 0 ],
         [ -w, n ]]:

#n := h;
#e := h;
#s := h;
#w := h;

#assume(n > 0);
#assume(e > 0);
#assume(s > 0);
#assume(w > 0);
#interface(showassumed=0);

#t := mtaylor(f(x,y), [x,y], 3);
#tt := unapply(t, [x,y]);

tt := (x,y) -> P0 + fx*x + fy*y + fxx*x^2/2 + fxy*x*y + fyy*y^2/2:

rels := {}:
for i from 1 to nops(dirs) do
    stencilpointname := "P" || (convert(i, string));
    rels := { op(rels), tt(dirs[i][1], dirs[i][2]) = parse(stencilpointname) };
end do:

# rels contains 8 equations for 5 unknowns.
#
# assume n,e,s,w are given
#
# find those fx, fy, fxx, fxy, fyy that most nearly fulfill all 8 equations.

# Solve
# A x = b,
#
# where x contains the five unknowns
# A contains the coefficients in front of the unknowns
# b contains the constants P0 ... P8

unknowns := [ fx, fy, fxx, fxy, fyy ]:
b := { op(unknowns) }:
#b := < P0, P1, P2, P3, P4, P5, P6, P7, P8 >;
#b := { P0, P1, P2, P3, P4, P5, P6, P7, P8 };

#solve(rels, unknowns);
with(LinearAlgebra):
res := LeastSquares(rels, b) assuming n > 0, e > 0, s > 0, w > 0:

for i from 1 to nops(res) do
    member(lhs(res[i]), unknowns, 'k'); # results may occur out-of-order
    for j from 0 to nops(dirs) do
        stencilpointname := "P" || (convert(j, string));
        c := simplify(coeff(rhs(res[i]), parse(stencilpointname)));
        code := CodeGeneration[C](c, output = string,
                                  deducetypes = false, resultname = "c");
        printf(" M0[%d+NDO*k+(%d+NSP*k)*(2*NDO)] = %s\n", (k-1), j, substring(code, 5..-1));
    end do:
end do:


