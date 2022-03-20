\\ Script to discard Forms 34, 35 and 36 in PARI/GP using Mazur's trick with q=5,7.

\\=============================================================================================================

\\ Form 34:

P=read("form34Pol.txt");
a5=read("form34Coef5.txt");
a7=read("form34Coef7.txt");
P=subst(P,x,a);

v5=[-4, -2, 0, 2, 4]  
w5=vector(length(v5),k,norm(-Mod(a5,P)-v5[k]));
w5=concat(w5,norm(Mod(a5,P)^2-36));
v7=[-14,-12,-6,4] 
w7=vector(length(v7),k,norm(Mod(a7,P)^2-v7[k]*(-1)-2*7*(-1)));
Primes=[];for(i=1,length(w5),for(j=1,length(w7),Primes=concat(Primes,factor(gcd(w5[i],w7[j]))[,1]~)))

Set(Primes)

\\=============================================================================================================

\\ Form 35:

P=read("form35Pol.txt");
a5=read("form35Coef5.txt");
a7=read("form35Coef7.txt");
P=subst(P,x,a);

v5=[-4, -2, 0, 2, 4]  
w5=vector(length(v5),k,norm(-Mod(a5,P)-v5[k]));
w5=concat(w5,norm(Mod(a5,P)^2-36));
v7=[-14,-12,-6,4] 
w7=vector(length(v7),k,norm(Mod(a7,P)^2-v7[k]*(-1)-2*7*(-1)));
Primes=[];for(i=1,length(w5),for(j=1,length(w7),Primes=concat(Primes,factor(gcd(w5[i],w7[j]))[,1]~)))

Set(Primes)

\\=============================================================================================================

\\ Form 36:

P=read("form36Pol.txt");
a5=read("form36Coef5.txt");
a7=read("form36Coef7.txt");
P=subst(P,x,a);

v5=[-4, -2, 0, 2, 4]  
w5=vector(length(v5),k,norm(-Mod(a5,P)-v5[k]));
w5=concat(w5,norm(Mod(a5,P)^2-36));
v7=[-14,-12,-6,4] 
w7=vector(length(v7),k,norm(Mod(a7,P)^2-v7[k]*(-1)-2*7*(-1)));
Primes=[];for(i=1,length(w5),for(j=1,length(w7),Primes=concat(Primes,factor(gcd(w5[i],w7[j]))[,1]~)))

Set(Primes)
