All the information compiled by F.Bello/E.Pierzchala on The fortran part of multiscat can be found in 
the HTML docs. This represents the structure of multiscat AFTER editing by F.Bello/E.Pierzchala.

Fay's guide is included as most of the information there is still relevant, although some is outdated.
Where Fay's guide contradicts the HTML docs consider the HTML to be correct as they are more up to date.

The structure of the HTML docs is not perfect and some the information is in awkward places: 
try clicking everything to get familiar with the structure, a lot of the names act like links 
to info pages. Under the "DependancyTree" there is an overall inheritance diagram as well as a 
description for each function and file. The "FlowChart" contains a diagram showing hoe the functions
are executed in sequency when the code is run. 

To edit the ULM documents the software "StarUML" was used (can be found for free here 
http://staruml.io/ at the time of writing). This allows .mdj files to be opened. The documentation
can be exported from this software as HTML and then read in an internet browser.
	
	
Papers used in this program:
	- 1:
		D.E. Manolopoulos, R.E. Wyatt,
		Quantum scattering via the log derivative version of the Kohn variational principle,
		Chemical Physics Letters,
		Volume 152, Issue 1, 1988, Pages 23-32
		
	- 2: 
		D.E. Manolopoulos, R.E. Wyatt,
		Iterative Solution in Quantum Scattering Theory. The log Derivative Kohn Approach
		J. Chem. Soc., Faraday Trans., 1990,86, 1641-1648 

	- 3:
		Yousef Saad
		Iterative Methods for Sparse Linear Systems
		Second edition, 2000
		
Gaussian-Legendre Quadrature is fairly well known, Gausian-Lobatto quadrature less so. It uses
	Lobatto polynomials (equivalent to the derivative of Legendre Polynomials) to determine 
        points and weights to use.
        <see GaussianQuadrature.pdf, bottom of page 6>(author unkown)

Lobatto Quadrature needs to find roots of Legendre polynomial's drivative. This is done via a 
        non-standard method:
	-Roots of the derivative are approximately half way between sucessive roots of the polynomial
	-To find the roots of the Legendre polynomial it uses Tricomi's approach
	-Uses Bonnet's recursion formula to evaluate the Legendre polynomials
	-Uses the differential equations that define Legendre polynomials to find derivatives
	-Uses Newton Raphson root finding to find the zeros of the derivative
	 (these are the x values needed for Lobatto quadrature rule)
	-Evaluates wheights to be used at each point <see GaussianQuadrature.pdf, bottom of page 6>

	
Lobatto-shape functions and their usefulness are explained in section 3 of [1]. 
	(Plot of them fig.1 is incorrect, correct version is shown in fig. 1 of [2] )

In section 4 of [1] it's explained how they're used to determine entries in matrices 
	used for determining scattering.
	
GMRES is well explained in chapter 6 of [3], and choice for preconditioner in
	section 3.2 of [2]





