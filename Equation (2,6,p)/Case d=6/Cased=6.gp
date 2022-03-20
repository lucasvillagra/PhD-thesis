\\ Discarding new[i] for i=52,...,58 in Pari/GP

\\ new[52]

P=read("Form52Pol.txt");
coef5=read("Form52Coef5.txt");
coef11=read("Form52Coef11.txt");
P=subst(P,x,a);
Rt=nfisincl(polcyclo(4,a),P)[1];

n5=norm(coef5*Mod(1,P))*norm(coef5*Mod(1,P)+3*Rt)*norm(coef5*Mod(1,P)+6*Rt);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

factor(gcd(n5,n11))

\\ new[53] 

P=read("Form53Pol.txt");
coef5=read("Form53Coef5.txt");
coef11=read("Form53Coef11.txt");
P=subst(P,x,a);
Rt=nfisincl(polcyclo(4,a),P)[1];

n5=norm(coef5*Mod(1,P))*norm(coef5*Mod(1,P)+3*Rt)*norm(coef5*Mod(1,P)+6*Rt);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

factor(gcd(n5,n11))

\\ new[54] 

P=read("Form54Pol.txt");
coef5=read("Form54Coef5.txt");
coef11=read("Form54Coef11.txt");
P=subst(P,x,a);
Rt=nfisincl(polcyclo(4,a),P)[1];

n5=norm(coef5*Mod(1,P))*norm(coef5*Mod(1,P)+3*Rt)*norm(coef5*Mod(1,P)+6*Rt);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

factor(gcd(n5,n11))

\\ new[55] 

P=read("Form55Pol.txt");
coef5=read("Form55Coef5.txt");
coef11=read("Form55Coef11.txt");
P=subst(P,x,a);
Rt=nfisincl(polcyclo(4,a),P)[1];

n5=norm(coef5*Mod(1,P))*norm(coef5*Mod(1,P)+3*Rt)*norm(coef5*Mod(1,P)+6*Rt);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

factor(gcd(n5,n11))

\\ new[56] 

P=read("Form56Pol.txt");
coef5=read("Form56Coef5.txt");
coef11=read("Form56Coef11.txt");
P=subst(P,x,a);
Rt=nfisincl(polcyclo(4,a),P)[1];

n5=norm(coef5*Mod(1,P))*norm(coef5*Mod(1,P)+3*Rt)*norm(coef5*Mod(1,P)+6*Rt);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

factor(gcd(n5,n11))

\\ new[57]

P=read("Form57Pol.txt");
coef5=read("Form57Coef5.txt");
coef11=read("Form57Coef11.txt");
P=subst(P,x,a);

n5=norm(coef5*Mod(1,P))*norm(coef5^2*Mod(1,P)+8)*norm(coef5^2*Mod(1,P)+36);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

factor(gcd(n5,n11))

\\ new[58] 

P=read("Form58Pol.txt");
coef5=read("Form58Coef5.txt");
coef11=read("Form58Coef11.txt");
P=subst(P,x,a);
Rt=nfisincl(polcyclo(4,a),P)[1];

n5=norm(coef5*Mod(1,P))*norm(coef5*Mod(1,P)+3*Rt)*norm(coef5*Mod(1,P)+6*Rt);
n11=norm(coef11*Mod(1,P))*norm(coef11^2*Mod(1,P)-9)*norm(coef11^2*Mod(1,P)-36)*norm(144-coef11^2*Mod(1,P));

gcd(n5,n11)

