/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

//#include <cmath>
#include "mixGCDSM_Flux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mixGCDSM_Flux, 0);
addToRunTimeSelectionTable(LESModel, mixGCDSM_Flux, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void mixGCDSM_Flux::updateSubGridScaleFields()
{
//    Info<< "calling mixGCDSM_Flux_updateSubGridScaleFields" << endl;
    muSgs_ = ck_*rho()*sqrt(k_)*delta();
    muSgs_.correctBoundaryConditions();

//    Info<< "min/max(mu) = " << min(mu()).value() << ", " << max(mu()).value() << endl;
//    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_ = mu()/Pr_+muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixGCDSM_Flux::mixGCDSM_Flux
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, rho, U, phi, thermoPhysicalModel, turbulenceModelName),
    GenEddyVisc(rho, U, phi, thermoPhysicalModel),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    ckp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ckp",
            coeffDict_,
            0.008
        )
    ),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.05
        )
    ),

    Pr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Pr",
            coeffDict_,
            0.7
        )
    ),

    cet_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cet",
            coeffDict_,
            1.0
        )
    )
{
    updateSubGridScaleFields();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
tmp<volSymmTensorField> mixGCDSM_Flux::B() const
{
//    Info<< "calling mixGCDSM_Flux_B" << endl;

    tmp<volTensorField> Du = fvc::grad(U());
    tmp<volTensorField> DuT = T(Du());
    tmp<volSymmTensorField> BB = symm(DuT() & Du());

//    tmp<volSymmTensorField> BB = scaleSimilarity::B();

    tmp<volScalarField> tB;
    tB = BB().component(symmTensor::XX)+BB().component(symmTensor::YY)+BB().component(symmTensor::ZZ);
    dimensionedScalar zero("zero", dimensionSet(0,0,-2,0,0,0,0),scalar(1.0e-30)); 
    tB = tB + zero;
    BB = 2.0*(k_/tB)*BB;

////    tmp<volSymmTensorField> SB = (Smagorinsky::B());
//    tmp<volSymmTensorField> SB = - ckp_*sqrt(k_)*delta()*twoSymm(fvc::grad(U()));
//    BB = BB+SB;
   

    if (ckp_.value() > 0.000001) {
    tmp<volSymmTensorField> SB = ckp_*sqrt(k_)*delta()*delta()*delta()*symm(fvc::grad(U()));
    tmp<volScalarField> LapSM;

    LapSM = fvc::laplacian(SB().component(symmTensor::XX));
    BB().component(symmTensor::XX) = BB().component(symmTensor::XX)+LapSM;

    LapSM = fvc::laplacian(SB().component(symmTensor::XY));
    BB().component(symmTensor::XY) = BB().component(symmTensor::XY)+LapSM;

    LapSM = fvc::laplacian(SB().component(symmTensor::XZ));
    BB().component(symmTensor::XZ) = BB().component(symmTensor::XZ)+LapSM;

    LapSM = fvc::laplacian(SB().component(symmTensor::YY));
    BB().component(symmTensor::YY) = BB().component(symmTensor::YY)+LapSM;

    LapSM = fvc::laplacian(SB().component(symmTensor::YZ));
    BB().component(symmTensor::YZ) = BB().component(symmTensor::YZ)+LapSM;

    LapSM = fvc::laplacian(SB().component(symmTensor::ZZ));
    BB().component(symmTensor::ZZ) = BB().component(symmTensor::ZZ)+LapSM;
    }

    return (BB);
//    return ((2.0/3.0)*I)*k() - (muSgs_/rho())*dev(twoSymm(fvc::grad(U())));
}

tmp<volScalarField> mixGCDSM_Flux::epsilon() const
{
//    Info<< "calling mixGCDSM_Flux_epsilon" << endl;
//    tmp<volSymmTensorField> D = symm(fvc::grad(U()));
   
//    return (B()&&D);
//    return (GenEddyVisc::B()&&D);
//    return (GenEddyVisc::epsilon());
//    return 2*muEff()/rho()*magSqr(symm(fvc::grad(U())));
    return (2*mu()/rho()*magSqr(symm(fvc::grad(U())))+ce_*k_*sqrt(k_)/delta()); //this is the total dissipation rate: resolved part+sgs part
//    return (ce_*k_*sqrt(k_)/delta());
}

tmp<volSymmTensorField> mixGCDSM_Flux::devRhoBeff() const
{
//    Info<< "calling mixGCDSM_Flux_devRhoBReff" << endl;
    return rho()*dev(B());
//    return -muEff()*dev(twoSymm(fvc::grad(U())));
}

tmp<fvVectorMatrix> mixGCDSM_Flux::divDevRhoBeff
(
    volVectorField& U
) const
{
//    Info<< "calling mixGCDSM_Flux_divDevRhoBReff" << endl;
//    return
//    (
//        fvc::div(rho()*B() + 0.05*muSgs_*fvc::grad(U))
//      + fvc::laplacian(0.95*muSgs_, U, "laplacian(muEff,U)")
//      - fvm::laplacian(muEff(), U)
//      - fvc::div(mu()*dev2(T(fvc::grad(U))))
//    );
    return
    (
      fvc::div(rho()*B())
      - fvm::laplacian(mu(), U)
//      fvc::div(rho()*B())
//      - fvm::laplacian(mu(), U) - fvc::div(mu()*dev2(T(fvc::grad(U))))
    );
}

// nonlinear model, Hao LU 2017/04/27
tmp<fvScalarMatrix> mixGCDSM_Flux::divDevRhoQeff
(
   volVectorField& U, 
   volScalarField& theta
) const
{
    tmp<volVectorField> Ds = fvc::grad(theta);
//    tmp<volVectorField> DS = fvc::grad(Myscalar());
//    const volVectorField& DS = fvc::grad(Myscalar);

    tmp<volTensorField> Du = fvc::grad(U);
    tmp<volTensorField> DuT = T(Du());
//    tmp<volSymmTensorField> BB = symm(DuT() & Du());
//    tmp<volSymmTensorField> SS = symm((Du() + DuT())/2);
//    tmp<volScalarField> tB;
//    tB = BB().component(symmTensor::XX)+BB().component(symmTensor::YY)+BB().component(symmTensor::ZZ);
////    dimensionedScalar zero("zero", dimensionSet(0,0,-2,0,0,0,0),scalar(1.0e-30)); 
////    dimensionedScalar zero("zero", (U.dimensions()/dimLength)*(U.dimensions()/dimLength),scalar(1.0e-30)); 
////    dimensionedScalar zero("zero", (U.dimensions()/dimLength)*(U.dimensions()/dimLength),SMALL); 
////    dimensionedScalar zero("zero", tB().dimensions(),SMALL); 
//    dimensionedScalar zero("zero", tB().dimensions(),scalar(1.0e-30)); 
//    tB = tB + zero;
//    BB = (1/tB)*BB;
//    tmp<volScalarField> P = -BB() && SS();
    

    tmp<volVectorField> BT = DuT() & Ds();
    tmp<volScalarField> MagBT = mag(BT());
//    dimensionedScalar zeroBT("zeroBT", MagBT().dimensions(),SMALL); 
    dimensionedScalar zeroBT("zeroBT", MagBT().dimensions(),scalar(1.0e-30)); 
    MagBT = MagBT + zeroBT;
    BT = (1/MagBT)*BT;
    tmp<volScalarField> PT = -BT() & Ds();
    
//    tmp<volScalarField> magQ = pos(P())*pos(PT())*2.828427*delta()*delta()/CeCet_*PT()*P();
    tmp<volScalarField> magQ = pos(PT())*delta()/cet_*PT()*sqrt(2*k_);
    BT = magQ*BT;
    
    return 
    (
        fvc::div(rho()*BT())
      - fvm::laplacian(mu()/Pr_, theta)
    );
//    return
//    (
//      - fvm::laplacian(alphaEff(), theta)
//    );
}


void mixGCDSM_Flux::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);

    volScalarField divU(fvc::div(phi()/fvc::interpolate(rho())));
//    volScalarField G(GName(), 2*muSgs_*(gradU && dev(symm(gradU))));
    volScalarField G(GName(), -rho()*(B() && (symm(gradU))));

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho(), k_)
      + fvm::div(phi(), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::SuSp(2.0/3.0*rho()*divU, k_)
      - fvm::Sp(ce_*rho()*sqrt(k_)/delta(), k_)
    );

    kEqn().relax();
    kEqn().solve();

    bound(k_, kMin_);

    updateSubGridScaleFields();
}


bool mixGCDSM_Flux::read()
{
    if (GenEddyVisc::read())
    {
        ckp_.readIfPresent(coeffDict());

        ck_.readIfPresent(coeffDict());

        Pr_.readIfPresent(coeffDict());

        cet_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
