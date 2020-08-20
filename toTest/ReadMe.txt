B.E. Pierzchala 08/2020

Two files here are implementations of tshapes in multiscat/scatsub.
Tshape0 is the one currecntly being used by the program and Tshape1 is same function 
rewritten and clarified.

What happened in the process is that if tshapes is implemented the way it should be
acoriding to [1], output is totally diferent to the implementation currently used 
by order of 10^3. This led to putting in place dummy variable 'hello', to account
for different method used previously, which I traced analytically.

So technically letting hello be what it is now, gives the same result, yet the two
differ probably because of roundoff errors ( 10^-7 differences for n = 500 ),
and if hello is set to 1, or deleted - it gives the results that [1] uses...

I am keen to believe that replacing in scatsub tshapes by implementation in Tshape1 would be
beneficial, as the code is written more clearly and one can easier trace it back to paper,
making it easier for future users.

What is left to do, as we ran out of time, is to run numerical experiments to see :

1. The difference of results of multiscat between using Tshape0 and Tshape1,
	I would argue that Tshape1 has less roundoff errors as it does less operations.
	( Tshape0 has integrated square roots into elements, then multiplies them and divides)
	
2. If the results using Tshape1 with hello = 1, gives results closer to theory we know
	or experiments we have. If not, it should be investigated why T matrix has been
	redefined between paper and code based on that paper.


1:
	D.E. Manolopoulos, R.E. Wyatt,
	Quantum scattering via the log derivative version of the Kohn variational principle,
	Chemical Physics Letters,
	Volume 152, Issue 1, 1988, Pages 23-32