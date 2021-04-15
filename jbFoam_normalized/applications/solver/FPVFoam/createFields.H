Info<< "Creating combustion model\n" << endl;
//Info<< "combustion sucess 1\n" << endl;
autoPtr<combustionModels::rhoCombustionModel> combustion
(
    combustionModels::rhoCombustionModel::New
    (
        mesh
    )
);
Info<< "combustion sucess 2\n" << endl;
combustion->correct();

Info<< "combustion sucess\n" << endl;

rhoReactionThermo& thermo = combustion->thermo();
 
Info<< "thermo sucess\n" << endl;

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

//Info << "rho is " <<endl;
//Info << rho << endl;



Info<< "rho sucess 1\n" << endl;

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

// Mixture Fraction
volScalarField& Z = thermo.Z();

// NO
//volScalarField& NO = thermo.NO();

// FPV progress variable
volScalarField& PVmajor = thermo.PVmajor();

// FPV progress variable source
volScalarField& OMGmajor = thermo.OMGmajor();

// FPV progress variable
volScalarField& PVNO = thermo.PVNO();

// FPV progress variable source
volScalarField& OMGNO = thermo.OMGNO();

// FPV progress variable
volScalarField& NOTransport = thermo.NOTransport();

// FPV progress variable source
volScalarField& OMGNOTransport = thermo.OMGNOTransport();



#include "compressibleCreatePhi.H"

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt("dpdt", fvc::ddt(p));

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

const volScalarField& psi = thermo.psi();
