#!/usr/local/bin/maple -q

# Gram-Schmidt process for generating orthogonal polynomials
#

if not assigned(n) then n := 10; end if:

phi0 := 1:
phi1 := x - 1/2:

memos := Array([seq(0, x=1..1000)]):

Bk := proc(k)
    local last_phi, top, bottom;

    last_phi := phi_k(k-1);
    top    := int(x * last_phi^2, x=0..1);
    bottom := int(last_phi^2, x=0..1);

    top / bottom;
end proc:

Ck := proc(k)
    local last_phi, last_phi2, top, bottom;

    last_phi  := phi_k(k-1);
    last_phi2 := phi_k(k-2);

    top    := int(x * last_phi * last_phi2, x=0..1);
    bottom := int(last_phi2^2, x=0..1);

    top / bottom;
end proc:

do_phi := proc(n)
    local res;

    if n = 0 then
        phi0;
    elif n = 1 then
        phi1;
    else
        res := (x - Bk(n)) * phi_k(n-1) - Ck(n) * phi_k(n-2);
        simplify(expand(res));
    end if;
end proc:

phi_k := proc(n)
    global memos;

    if memos[n+1] = 0 then
        memos[n+1] := do_phi(n);
    end if;

    memos[n+1];
end proc:

genphis := proc(n)
    local i;
    for i from 1 to n do
        phi_k(i);
    end do;
end proc;

print(n);
genphis(n):

failures := 0:
for i from 1 to numelems(memos)-1 do
    for j from 1 to 9 do
        if i <> j then
            res := int(memos[i] * memos[i+1], x=0..1);
            if res <> 0 then
                printf("error. found phi %i and %j whos inner product is not zero\n", i, j);
                failures = failures + 1;
            end if;
        end if;
    end do;
end do:

if failures = 0 then
    printf("success: the set of polynomials forms an orthonormal space\n");
end if:

memos(4);
