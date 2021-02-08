
surfaceScalarField phiC =
(
	(
		( 1/H * (D1-D2)/(fvc::interpolate(alpha1) + fvc::interpolate(alpha2)/H) - (D1 - D2/H) )/(fvc::interpolate(alpha1) + fvc::interpolate(alpha2)/H )	
	) * fvc::snGrad(alpha1) 

) * mesh.magSf();


solve //fvScalarMatrix CEqn
 (
		  fvm::ddt(C) 
		+ fvm::div(phi,C)
		- fvm::laplacian(fvc::interpolate(D),C)
		- fvm::laplacian ( alpha1 * (D1-D2) * ( 1/(alpha1 + alpha2/H) -1 ) , C)
		- fvm::div(phiC,C,"div(phi,C)")
 );


