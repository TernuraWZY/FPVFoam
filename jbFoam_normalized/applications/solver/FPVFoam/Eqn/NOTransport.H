fvScalarMatrix NOTransportEqn
( 
    (
    fvm::ddt(rho, NOTransport)
    + fvm::div(phi, NOTransport)
    - fvm::laplacian(turbulence->DPVEff(), NOTransport)
    == rho*OMGNOTransport                          //进程变量反应源项，在YSLFmodel中赋值。DPVEff仿造DZEff
    )
);

NOTransportEqn.relax();
NOTransportEqn.solve(mesh.solver("NOTransport"));