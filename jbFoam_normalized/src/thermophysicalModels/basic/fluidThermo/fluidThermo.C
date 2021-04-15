/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/

#include "fluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::fluidThermo(const fvMesh& mesh, const word& phaseName)
:
    basicThermo(mesh, phaseName),
    
    Z_
    (
        IOobject
        (
            "Z",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    varZ_
    (
        IOobject
        (
            "varZ",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Chi_
    (
        IOobject
        (
            phasePropertyName("chi"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
	
    PVmajor_
    (
        IOobject
        (
            "PVmajor",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),  //add by Jiang 添加进程变量
    
        PVNO_
    (
        IOobject
        (
            "PVNO",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),  //add by Jiang 添加进程变量

    OMGmajor_
    (
        IOobject
        (
            "OMGmajor",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),     //add by Jiang 添加进程变量源项
    
    OMGNO_
    (
        IOobject
        (
            "OMGNO",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),     //add by Jiang 添加进程变量源项
    
    NO_
    (
        IOobject
        (
            "NO",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),     //add by wzy NO fields
        NOTransport_
    (
        IOobject
        (
            "NOTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),     //add by wzy NO transport
        OMGNOTransport_
    (
        IOobject
        (
            "OMGNOTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),    //add by wzy NO transport
        OMGNOTransportPos_
    (
        IOobject
        (
            "OMGNOTransportPos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),     //add by wzy NO transport
        OMGNOTransportNeg_
    (
        IOobject
        (
            "OMGNOTransportNeg",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    )     //add by wzy NO transport    
    

{}



Foam::fluidThermo::fluidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    basicThermo(mesh, dict, phaseName),

    
    Z_
    (
        IOobject
        (
            phasePropertyName("Z"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    varZ_
    (
        IOobject
        (
            phasePropertyName("varZ"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    Chi_
    (
        IOobject
        (
            phasePropertyName("chi"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    PVmajor_
    (
        IOobject
        (
            "PVmajor",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),  //add by Jiang 添加进程变量
    
        PVNO_
    (
        IOobject
        (
            "PVNO",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),  //add by Jiang 添加进程变量

    OMGmajor_
    (
        IOobject
        (
            "OMGmajor",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),    //add by Jiang 添加进程变量源项
    
    OMGNO_
    (
        IOobject
        (
            "OMGNO",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),    //add by Jiang 添加进程变量源项

    NO_
    (
        IOobject
        (
            "NO",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),    //add by wzy NO fields
            NOTransport_
    (
        IOobject
        (
            "NOTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),     //add by wzy NO transport
        OMGNOTransport_
    (
        IOobject
        (
            "OMGNOTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),     //add by wzy NO transport
        OMGNOTransportPos_
    (
        IOobject
        (
            "OMGNOTransportPos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    ),     //add by wzy NO transport
        OMGNOTransportNeg_
    (
        IOobject
        (
            "OMGNOTransportNeg",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, -1, 0, 0, 0, 0)
    )     //add by wzy NO transport
    
    

{}     //初步认为通过字典的数据交换的各个变量是.H中的常量，不能变化，只是用来查表，而网格点直接数据交换的变量是非常量，用来迭代。


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}


Foam::volScalarField& Foam::fluidThermo::Z()
{
    return Z_;
}

const Foam::volScalarField& Foam::fluidThermo::Z() const
{
    return Z_;
}

Foam::volScalarField& Foam::fluidThermo::varZ()
{
    return varZ_;
}


const Foam::volScalarField& Foam::fluidThermo::varZ() const
{
    return varZ_;
}

Foam::volScalarField& Foam::fluidThermo::Chi()
{
    return Chi_;
}

const Foam::volScalarField& Foam::fluidThermo::Chi() const
{
    return Chi_;
}

//  Non-const access allowed for transport equations
Foam::volScalarField& Foam::fluidThermo::PVmajor()
{
    return PVmajor_;
}

//- progress valiable
const Foam::volScalarField& Foam::fluidThermo::PVmajor() const
{
    return PVmajor_;
}              //add Jiang 给进程变量定义两种类型，分别用来插值和迭代

Foam::volScalarField& Foam::fluidThermo::PVNO()
{
    return PVNO_;
}

//- progress valiable
const Foam::volScalarField& Foam::fluidThermo::PVNO() const
{
    return PVNO_;
}              //add Jiang 给进程变量定义两种类型，分别用来插值和迭代

//  Non-const access allowed for transport equations
Foam::volScalarField& Foam::fluidThermo::OMGmajor()
{
    return OMGmajor_;
}

//- progress valiable

Foam::volScalarField& Foam::fluidThermo::OMGNO()
{
    return OMGNO_;
}

Foam::volScalarField& Foam::fluidThermo::NO_Field()
{
    return NO_;
}

Foam::volScalarField& Foam::fluidThermo::NOTransport()
{
    return NOTransport_;
}

Foam::volScalarField& Foam::fluidThermo::OMGNOTransport()
{
    return OMGNOTransport_;
}

Foam::volScalarField& Foam::fluidThermo::OMGNOTransportPos()
{
    return OMGNOTransportPos_;
}

Foam::volScalarField& Foam::fluidThermo::OMGNOTransportNeg()
{
    return OMGNOTransportNeg_;
}

// ************************************************************************* //
