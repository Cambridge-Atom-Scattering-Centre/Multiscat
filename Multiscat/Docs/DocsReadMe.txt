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

Uses Lobatto rule to approximate integral ()

Lobatto needs to find roots of Legendre polynomial (we could use the monic Legendre polynomials):
	-Uses Bonnet's recursion formula to evaluate the polynomials, 
	see RecursiveLegendre.pdf (source: University of souther california CACS)
	-To find the roots it uses Tricomi's approach 
	-Followed by Newton Raphson ()