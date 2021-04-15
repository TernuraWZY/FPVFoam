/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "WALE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WALE, 0);
addToRunTimeSelectionTable(LESModel, WALE, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void WALE::updateSubGridScaleFields(const volTensorField& gradU)
{
    volSymmTensorField Sd(dev(symm(gradU & gradU)));
    volScalarField k = sqr(sqr(cw_)*delta()/ck_)*
            (
                pow3(magSqr(Sd))
               /(
                   sqr
                   (
                       pow(magSqr(symm(gradU)), 5.0/2.0)
                     + pow(magSqr(Sd), 5.0/4.0)
                   )
                 + dimensionedScalar
                   (
                       "SMALL",
                       dimensionSet(0, 0, -10, 0, 0),
                       SMALL
                   )
               )
            );
   muSgs_ = ck_*rho()*sqrt(k)*delta();
   muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WALE::WALE
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


    ck_			
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),

    cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cw",
            coeffDict_,
            0.325
        )
    )
{
    updateSubGridScaleFields(fvc::grad(U));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WALE::correct(const tmp<volTensorField>& gradU)
{

    GenEddyVisc::correct(gradU);
								
    updateSubGridScaleFields(gradU());	

}


bool WALE::read()
{
    if (GenEddyVisc::read())
    {
        ck_.readIfPresent(coeffDict());

        cw_.readIfPresent(coeffDict());

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
