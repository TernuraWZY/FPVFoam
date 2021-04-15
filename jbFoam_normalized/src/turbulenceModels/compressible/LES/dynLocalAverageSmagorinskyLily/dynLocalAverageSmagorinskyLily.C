/*---------------------------------------------------------------------------*\
dynLocalAverageSmagorinsky - Implementation of the dynamic Smagorinsky
    			     SGS model.
    
Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    
License
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dynLocalAverageSmagorinskyLily.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynLocalAverageSmagorinsky, 0);
addToRunTimeSelectionTable(LESModel, dynLocalAverageSmagorinsky, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dynLocalAverageSmagorinsky::updateSubGridScaleFields(const volTensorField& gradU, volScalarField& LijMij_, volScalarField& MklMkl_)
{
    Info << "dyn success " <<endl;
    volSymmTensorField D(symm(gradU));

    // The SGS viscosity is bounded so that muEff cannot become negative.
    // Values are limited here, and not in muEff, for consistency in stored
    // data and in submodels using muSgs().
    // No warning message is printed when this limitation is applied.
    muSgs_ = max(rho()*cD_(gradU, LijMij_, MklMkl_)*sqr(delta())*sqrt(2.0*(D && D)), -mu());
    muSgs_.correctBoundaryConditions();
    //alphaZmixSgs_ =alphaZmix();
    //alphaZmixSgs_.correctBoundaryConditions();
    //alphaProgSgs_ = alphaProg();
    //alphaProgSgs_.correctBoundaryConditions();
    //ZvarCoeff_ = max(cDzv_(), 0.0);
    //ZvarCoeff_.correctBoundaryConditions();
    Info << "\nDyn update success\n" << endl;
}
 

volScalarField dynLocalAverageSmagorinsky::cD_(const volTensorField& gradU, volScalarField& LijMij_, volScalarField& MklMkl_) const
{
    volSymmTensorField D(symm(gradU));
    volScalarField magSij = sqrt(2.0*(D && D));
    
    volSymmTensorField Mij = 2.0*sqr(delta())*(filter_(D*rho()*magSij)-4.0*filter_(rho())*filter_(magSij)*filter_(D));
    //volSymmTensorField Mij = 2.0*sqr(delta())*(filter_(dev(D)*rho()*magSij)-4.0*filter_(rho())*filter_(magSij)*(filter_(D)-filter_(tr(D))*I/3));
    //volSymmTensorField Mij = 2*sqr(delta())*(filter_(D*rho()*magSij)-2*filter_(rho())*filter_(magSij)*filter_(D));
 
 //   volSymmTensorField Bij = -2*sqr(2*delta())*filter_(rho())*filter_(magSij)*(filter_(D)-filter_(tr(D))*I/3);
 //   volSymmTensorField Aij = -2*sqr(delta())*rho()*magSij*dev(D);

    volSymmTensorField Lij = filter_(rho()*sqr(U())) - sqr(filter_(rho()*U()))/filter_(rho());
    //volSymmTensorField Mij = Bij - filter_(Aij);

    volScalarField LijMij = Lij && Mij;
    volScalarField MklMkl = Mij && Mij;

    LijMij_ = LijMij;
    MklMkl_ = MklMkl;
    
    // Performing local average on cell faces on return
    // return fvc::average(LLMM)/MMMM;
    return fvc::average(LijMij)/fvc::average(MklMkl);
}


/*volScalarField dynLocalAverageSmagorinsky::cDscZmix_(const volTensorField& gradU) const
{
    volSymmTensorField D(symm(gradU));
    volScalarField magSij = sqrt(2*(D && D));
    //volVectorField Mk = 2.0*sqr(delta())*(filter_(rho()*magSij*fvc::grad(Zmix()))-
    //                            2.0*filter_(rho())*filter_(magSij)*filter_(fvc::grad(Zmix()))); 
    volVectorField Mk = 2.0*sqr(delta())*(filter_(rho()*magSij*fvc::grad(Zmix()))-
                               4.0*filter_(rho())*filter_(magSij)*filter_(fvc::grad(Zmix()))); 
    volVectorField Li = filter_(rho()*U()*Zmix()) - filter_(rho()*U())*filter_(rho()*Zmix())/filter_(rho());
    volScalarField MkMk = Mk & Mk;
    volScalarField LiMk = Li & Mk;
    
    
    return (fvc::average(LiMk))/(fvc::average(MkMk));

//    volVectorField Li = filter_(rho()*U()*gradZmix()) - (filter_(rho())*filter_(U)*filter_(Zmix));
//    volVectorField Mi = filter_(rho()*magSij*gradZmix()) - ratio*ratio*filter_(rho)*filter_(magSij)*filter_(gradZmix());

//    volScalarField LiMi = Li & Mi;
//    volScalarField MkMk = Mi & Mi;

 //   LiMi_ = LiMi;
 //   MkMk_ = MkMk;
    
    // Performing local average on cell faces on return
    // return fvc::average(LLMM)/MMMM;
 //   return fvc::average(LiMi)/fvc::average(MkMk);
}*/


/*volScalarField dynLocalAverageSmagorinsky::cDscProg_(const volTensorField& gradU) const
{
    volSymmTensorField D(symm(gradU));
    volScalarField magSij = sqrt(2*(D && D));
    //volVectorField Mk = 2.0*sqr(delta())*(filter_(rho()*magSij*fvc::grad(Zmix()))-
    //                            2.0*filter_(rho())*filter_(magSij)*filter_(fvc::grad(Zmix()))); 
    volVectorField Mk = 2.0*sqr(delta())*(filter_(rho()*magSij*fvc::grad(Prog()))-
                               4.0*filter_(rho())*filter_(magSij)*filter_(fvc::grad(Prog()))); 
    volVectorField Li = filter_(rho()*U()*Prog()) - filter_(rho()*U())*filter_(rho()*Prog())/filter_(rho());
    volScalarField MkMk = Mk & Mk;
    volScalarField LiMk = Li & Mk;
    
    
    return (fvc::average(LiMk))/(fvc::average(MkMk));

//    volVectorField Li = filter_(rho()*U()*gradZmix()) - (filter_(rho())*filter_(U)*filter_(Zmix));
//    volVectorField Mi = filter_(rho()*magSij*gradZmix()) - ratio*ratio*filter_(rho)*filter_(magSij)*filter_(gradZmix());

//    volScalarField LiMi = Li & Mi;
//    volScalarField MkMk = Mi & Mi;

 //   LiMi_ = LiMi;
 //   MkMk_ = MkMk;
    
    // Performing local average on cell faces on return
    // return fvc::average(LLMM)/MMMM;
 //   return fvc::average(LiMi)/fvc::average(MkMk);
}*/



/*volScalarField dynLocalAverageSmagorinsky::cDzv_(void) const
{
    volScalarField magGZ = mag(fvc::grad(Zmix()) & fvc::grad(Zmix()));
    volScalarField Mk = sqr(delta())*4.0*filter_(rho())*filter_(magGZ);
    volScalarField Li = filter_(rho()*Zmix()*Zmix()) - filter_(rho())*sqr(filter_(rho()*Zmix())/filter_(rho()));
    volScalarField MkMk = Mk * Mk;
    volScalarField LiMk = Li * Mk;


    return (fvc::average(LiMk))/(fvc::average(MkMk));

}*/


volScalarField dynLocalAverageSmagorinsky::cI_(const volTensorField& gradU) const
{
    volSymmTensorField D(symm(gradU));

    //volScalarField KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));
    volSymmTensorField Lij = filter_(rho()*sqr(U())) - sqr(filter_(rho()*U()))/filter_(rho());
    volScalarField Lkk = tr(Lij); 
    volScalarField magSij = sqrt(2.0*(D && D));
    volScalarField beta = 2.0*sqr(2*delta())*filter_(rho())*sqr(filter_(magSij));
    volScalarField alpha = filter_(2.0*rho()*sqr(delta())*sqr(magSij));

    // volScalarField mm =
    // sqr(delta())*(4*sqr(mag(filter_(D))) - filter_(sqr(mag(D))));

    // Locally averaging mmmm on cell faces
    // volScalarField mmmm = fvc::average(magSqr(mm));
    // mmmm.max(VSMALL);

    // Performing local average on cell faces on return
    return (fvc::average(Lkk))/(fvc::average(beta - alpha));
}

/*volScalarField dynLocalAverageSmagorinsky::cQ_(const volSymmTensorField& D) const
{
    const volScalarField& T = db().lookupObject<volScalarField>("T");
    volScalarField magSij = sqrt(2*(D && D));
    volVectorField Tk = -2*sqr(2*delta())*filter_(rho())*filter_(magSij)*fvc::grad(filter_(T)) +
                         sqr(delta())*filter_(rho()*magSij*fvc::grad(T));
    volVectorField Kj = filter_(rho()*U()*T) - filter_(rho()*U())*filter_(rho()*T)/filter_(rho());
    volScalarField TkTk = Tk & Tk;
    volScalarField KjTj = Kj & Tk;
    return (fvc::average(TkTk))/(fvc::average(KjTj));
}*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynLocalAverageSmagorinsky::dynLocalAverageSmagorinsky
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

    LijMij_
    (
        IOobject
        (
            "LijMij",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("LijMij", dimensionSet(2,-2,-4,0,0,0,0), 0.0)
    ),


    MklMkl_
    (
        IOobject
        (
            "MklMkl",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("MklMkl", dimensionSet(2,-2,-4,0,0,0,0), 0.0)
    ),
    

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    updateSubGridScaleFields(fvc::grad(U), LijMij_, MklMkl_);

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynLocalAverageSmagorinsky::correct(const tmp<volTensorField>& tgradU)
{
   // LESModel::correct(tgradU);
    const volTensorField& gradU = tgradU();
    GenEddyVisc::correct(gradU);
    volSymmTensorField D = symm(gradU);
    volScalarField magSij = sqrt(2.0*(D && D));
    k_ = 2*cI_(gradU)*sqr(delta()*magSij);
    updateSubGridScaleFields(gradU, LijMij_, MklMkl_);
}

bool dynLocalAverageSmagorinsky::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

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
