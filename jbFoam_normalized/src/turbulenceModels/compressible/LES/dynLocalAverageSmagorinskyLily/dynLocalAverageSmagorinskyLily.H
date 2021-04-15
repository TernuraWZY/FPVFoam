/*---------------------------------------------------------------------------*\
dynLocalAverageSmagorinsky - Implementation of the dynamic Smagorinsky
			     SGS model.
				 
Copyright Information
    Copyright (C) 1991-2009 OpenCFD Ltd.
    Copyright (C) 2010 Alberto Passalacqua 

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

Class
    Foam::compressible::LESModels::dynLocalAverageSmagorinsky

Description
    The dynamic Smagorinsky model for compressible flows.
     From Subgrid-Scale Models for Compressible LES,Martin, Piomelli, & Candler
     Theoretical Comput. Fluid Dynamics 2000 13:361-376


SourceFiles
    dynLocalAverageSmagorinsky.C
    
Author
    Ayodeji Ojofeitimi, ayodejio@umd.edu

Notes
    Implementation of the dynamic Smagorinsky model with coefficients cD and
    cI computed as local average of their face values to avoid numerical 
    instabilities. 

    Negative values of the effective viscosity are removed by clipping it to
    zero (muSgs is clipped to -mu)

    The code is known to work with OpenFOAM 1.7.x

\*---------------------------------------------------------------------------*/

#ifndef dynLocalAverageSmagorinsky_H
#define dynLocalAverageSmagorinsky_H

#include "GenEddyVisc.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynLocalAverageSmagorinsky Declaration
\*---------------------------------------------------------------------------*/

class dynLocalAverageSmagorinsky
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;
        volScalarField LijMij_;
        volScalarField MklMkl_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volTensorField& gradU, volScalarField& LijMij_, volScalarField& MklMkl_);
	
        //- Calculate coefficients cD, cI from filtering velocity field
        volScalarField cD_(const volTensorField& gradU, volScalarField& LijMij_, volScalarField& MklMkl_) const;
        volScalarField cI_(const volTensorField& gradU) const;
        //volScalarField cDsc_(const volTensorField& gradU) const;
        //volScalarField cDscZmix_(const volTensorField& gradU) const;
        //volScalarField cDscProg_(const volTensorField& gradU) const;
        // volScalarField cDzv_(void) const;
        // volScalarField cQ_(const volTensorField& gradU) const;
        

        // Disallow default bitwise copy construct and assignment
        dynLocalAverageSmagorinsky(const dynLocalAverageSmagorinsky&);
        dynLocalAverageSmagorinsky& operator=(const dynLocalAverageSmagorinsky&);


public:

    //- Runtime type information
    TypeName("dynLocalAverageSmagorinskyLily");

    // Constructors

        //- Construct from components
        dynLocalAverageSmagorinsky
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

    //- Destructor
    virtual ~dynLocalAverageSmagorinsky()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
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
