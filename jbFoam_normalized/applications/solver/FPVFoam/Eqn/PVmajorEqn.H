fvScalarMatrix PVmajorEqn
(
  (
      fvm::ddt(rho, PVmajor)
    + fvm::div(phi, PVmajor)
    - fvm::laplacian(turbulence->DPVEff(), PVmajor)
    == rho*OMGmajor                           //进程变量反应源项，在YSLFmodel中赋值。DPVEff仿造DZEff
  )
);

PVmajorEqn.relax();
PVmajorEqn.solve(mesh.solver("PVmajor"));