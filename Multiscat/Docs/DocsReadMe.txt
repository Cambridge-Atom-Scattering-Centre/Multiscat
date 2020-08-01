All the information compiled by F.Bello/E.Pierzchala on The fortran part of multiscat can be found in 
the HTML docs. This represents the structure of multiscat BEFORE editing by F.Bello/E.Pierzchala.

The structure of the HTML docs is not perfect and some the information is in awkward places: 
try clicking everything to get familiar with the structure, a lot of the names act like links 
to info pages

To edit the ULM documents the software "StarUML" was used (can be found for free here 
http://staruml.io/ at the time of writing). This allows .mdj files to be opened. The documentation
can be read in an internet browser once exported from this siftware as HTML.

Solver:
(where does gmres fit?)
(what are we integrating?)

Uses Gauss-Lobatto quadrature to approximate the integral 
<see GaussianQuadrature.pdf, bottom of page 6>(author unkown)

Lobatto needs to find roots of Legendre polynomial's drivative:
	-Roots of the derivative are approximately half way between sucessive roots of the polynomial
	-To find the roots it uses Tricomi's approach
	-Uses Bonnet's recursion formula to evaluate the Legendre polynomials
	-Uses the differential equations that define Legendre polynomials to find derivatives
	-Uses Newton Raphson root finding the zeros of the derivative
	 (these are the x values needed for Lobatto quadrature rule)
	-Evaluates wheights to be used at each point <see GaussianQuadrature.pdf, bottom of page 6>