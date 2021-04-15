fvScalarMatrix PVNOEqn
(
	(
		  fvm::ddt(rho, PVNO)
		+ fvm::div(phi, PVNO)
        - fvm::laplacian(turbulence->DPVEff(), PVNO)
		== rho*OMGNO
    )
);

PVNOEqn.relax();
PVNOEqn.solve(mesh.solver("PVNO"));