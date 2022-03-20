load "Mazur26p.mg";

/* Character chi */

Chi:= function(p)
if KroneckerSymbol(-6,p) eq -1 then
	f:=order(KroneckerSymbol(-1,p)*KroneckerSymbol(p,2)*KroneckerSymbol(p,3));
else
   f:=2*order(KroneckerSymbol(-1,p)*KroneckerSymbol(p,3));
end if;
return f;
end function;

/* The space */

G:=DirichletGroup(2^8*3^5);
eps:=Generators(G)[1]*Generators(G)[3];
M:=ModularSymbols(eps,2,1);
S:=NewSubspace(CuspidalSubspace(M));
new:=NewformDecomposition(S);
CM:=FormsWithCM(new);

/* Mazur's trick for forms without CM until the place 51 */

for i in [1..51] do
if i notin CM then
print(DiscardPlace(6,eps,Chi,new,i,1,13));  
end if;           
end for;

/* To discard forms new[i] for i=52,...,28 */

ApCand(6,5);
ApCand(6,11);

/* new[52] */

f:=PowerSeries(new[52],12);
Write("Form52Coef5.txt",Coefficient(f,5));
Write("Form52Coef11.txt",Coefficient(f,11));
Write("Form52Pol.txt",Parent(Coefficient(f,5)));

/* new[53] */

f:=PowerSeries(new[53],12);
Write("Form53Coef5.txt",Coefficient(f,5));
Write("Form53Coef11.txt",Coefficient(f,11));
Write("Form53Pol.txt",Parent(Coefficient(f,5)));

/* new[54] */

f:=PowerSeries(new[54],12);
Write("Form54Coef5.txt",Coefficient(f,5));
Write("Form54Coef11.txt",Coefficient(f,11));
Write("Form54Pol.txt",Parent(Coefficient(f,5)));

/* new[55] */

f:=PowerSeries(new[55],12);
Write("Form55Coef5.txt",Coefficient(f,5));
Write("Form55Coef11.txt",Coefficient(f,11));
Write("Form55Pol.txt",Parent(Coefficient(f,5)));

/* new[56] */

f:=PowerSeries(new[56],12);
Write("Form56Coef5.txt",Coefficient(f,5));
Write("Form56Coef11.txt",Coefficient(f,11));
Write("Form56Pol.txt",Parent(Coefficient(f,5)));

/* new[57] */

f:=PowerSeries(new[57],12);
Write("Form57Coef5.txt",Coefficient(f,5));
Write("Form57Coef11.txt",Coefficient(f,11));
Write("Form57Pol.txt",Parent(Coefficient(f,5)));

/* new[58] */

f:=PowerSeries(new[58],12);
Write("Form58Coef5.txt",Coefficient(f,5));
Write("Form58Coef11.txt",Coefficient(f,11));
Write("Form58Pol.txt",Parent(Coefficient(f,5)));
