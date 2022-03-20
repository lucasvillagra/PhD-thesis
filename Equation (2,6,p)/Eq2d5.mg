load "Mazur26p.mg";


/* Character chi */

Chi:= function(p)
Eps1:=Generators(DirichletGroup(15,CyclotomicField(4)));
Eps1:=Eps1[1]*Eps1[2];
Eps2:=Generators(DirichletGroup(20,CyclotomicField(4)));
Eps2:=Eps2[1]*Eps2[2];
if KroneckerSymbol(-5,p) eq -1 then
   f:=order(Eps1(p));
else
	f:=order(Eps2(p))*2;
end if;
return f;
end function;

/* First space */

G:=DirichletGroup(2^4*3^2*5^2,CyclotomicField(4));                           
eps:=(Elements(G)[50]);
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

/* Mazur's trick for forms without CM */

for i in [1..#new] do         
if i notin CM then
print(DiscardPlace(5,eps,Chi,new,i,13,50));
end if;
end for; 

/* Second space */

eps:=Extend(eps,2^4*3^3*5^2);
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

/* Mazur's trick for forms without CM */

for i in [1..#new] do         
if i notin CM then
print(DiscardPlace(5,eps,Chi,new,i,13,50));
end if;
end for; 

