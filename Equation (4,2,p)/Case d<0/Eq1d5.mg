load "Mazur42p.mg";

/* Character chi */

Chi:= function(p)
G:=DirichletGroup(2^7*5^2,CyclotomicField(4));                       
eps:=(Elements(G)[26]); 
f:=order(eps(p)*KroneckerSymbol(p,5)*KroneckerSymbol(p,2));
if KroneckerSymbol(-5,p) eq 1 then
   f:=2*order(eps(p));
end if;
return f;
end function;

/* First space */

print("Forms in Space 2^7*5^2:");
G:=DirichletGroup(2^7*5^2,CyclotomicField(4));                       
eps:=(Elements(G)[26]); 
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

print "There are", #new, "forms,", #CM, "of them having complex multiplication.";
print "Forms with CM:", CM;

print("Forms with CM:");
print(CM);

/* Mazur's trick for forms without CM */

print("Primes obtained via Mazur's trick for non-CM forms:");
for i in [1..#new] do
if i notin CM then
print(DiscardPlace(5,eps,Chi,new,i,1,30));
end if;
end for;

/* Second space */

print("Forms in Space 2^8*5^2:");
G:=DirichletGroup(2^8*5^2,CyclotomicField(4));                       
eps:=(Elements(G)[26]);   
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

print "There are", #new, "forms,", #CM, "of them having complex multiplication.";
print "Forms with CM:", CM;

print("Forms with CM:");
print(CM);

/* Mazur's trick for forms without CM */

print("Primes obtained via Mazur's trick for non-CM forms:");
for i in [1..#new] do
if i notin CM then
print(DiscardPlace(5,eps,Chi,new,i,1,20));
end if;
end for;

/* Mazur's trick for forms with CM */

print("Primes obtained via Mazur's trick for CM forms:");
for i in CM do                                 
print(MazurTrickMultiplicative(new,3,eps,[i]));
end for;
