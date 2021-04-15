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

#ifndef compressibleWALE_H
#define compressibleWALE_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class WALE Declaration
\*---------------------------------------------------------------------------*/

class WALE
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;     

        dimensionedScalar ck_;

	dimensionedScalar cw_;

//	volScalarField sMut;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volTensorField& gradU);

        // Disallow default bitwise copy construct and assignment
        WALE(const WALE&);
        WALE& operator=(const WALE&);


public:

    //- Runtime type information
    TypeName("WALE");


    // Constructors

        //- Constructor from components
        WALE
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName

        );


    // Destructor
    virtual ~WALE()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        //  calculated from the given velocity gradient
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const
	{
	  volSymmTensorField Sd(dev(symm(gradU & gradU)));
	  return sqr(sqr(cw_)*delta()/ck_)*
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
        }
        
        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(U()));
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
