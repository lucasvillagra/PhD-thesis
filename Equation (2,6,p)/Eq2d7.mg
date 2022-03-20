load "Mazur26p.mg";

/* Character chi */

Chi:= function(p)
if KroneckerSymbol(-7,p) eq -1 then
   f:=1;
else
	f:=2*order(KroneckerSymbol(p,3)*KroneckerSymbol(p,7));
end if;
return f;
end function;

/* First space */

G:=DirichletGroup(2*3*7^2);
eps:=Generators(G)[1]*Generators(G)[2];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);

/* Mazur's trick */

for i in [1..#new] do         
print(DiscardPlace(7,eps,Chi,new,i,13,23));
end for;      


/* Second space */

G:=DirichletGroup(2^2*3*7^2);
eps:=Generators(G)[2]*Generators(G)[3];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

/* Mazur's trick for forms without CM*/

for i in [2..#new] do
if i notin CM then
print(DiscardPlace(7,eps,Chi,new,i,1,20));
end if;
end for;

/* Third space */

G:=DirichletGroup(2*3^3*7^2);
eps:=Generators(G)[1]*Generators(G)[2];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);

/* Mazur's trick */

for i in [1..#new] do
DiscardPlace(7,eps,Chi,new,i,1,13);
end for;

/* Fourth space */

G:=DirichletGroup(2^2*3^3*7^2);
eps:=Generators(G)[2]*Generators(G)[3];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

/* Mazur's trick for forms without CM */

for i in [1..#new] do
if i notin CM then
print(DiscardPlace(7,eps,Chi,new,i,1,20));  
end if;           
end for;

