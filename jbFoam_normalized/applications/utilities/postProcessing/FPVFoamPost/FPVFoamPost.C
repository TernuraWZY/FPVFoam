/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Calculates species mass fractions and thermodynamic properties
    from given Z, chi and Zeta fields

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "tableSolver.H"

#include "fvCFD.H"
#include "rhoCombustionModel.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be reconstructed. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    Info<< "Creating combustion model\n" << endl;

    autoPtr<combustionModels::rhoCombustionModel> combustion
    (
        combustionModels::rhoCombustionModel::New
        (
            mesh
        )
    );

    rhoReactionThermo& thermo = combustion->thermo();    

    const IOdictionary combProps
    (
        IOobject
        (
            "combustionProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const IOdictionary tableProps
    (
        IOobject
        (
            "tableProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word combTypeName = combProps.lookup("combustionModel");
    const label tempOpen = combTypeName.find('<');
    const word modelType = combTypeName(0, tempOpen);
    dictionary coeffs_(combProps.subDict(modelType + "Coeffs"));

    Switch useScalarDissipation_(coeffs_.lookup("useScalarDissipation"));
    Switch useMixtureFractionVariance_(coeffs_.lookup("useMixtureFractionVariance"));

    scalarList chi_Param(tableProps.lookup("chi_param"));
    //scalar chiLimiter = max(chi_Param);

    scalarList x(3);
    scalarList y(3);
    List<List<int> > ubMajorIF(mesh.cells().size());
    List<scalarList> posMajorIF(mesh.cells().size());
    List<List<int> > ubMajorP(mesh.faces().size());
    List<scalarList> posMajorP(mesh.faces().size());

    List<List<int> > ubNOIF(mesh.cells().size());
    List<scalarList> posNOIF(mesh.cells().size());
    List<List<int> > ubNOP(mesh.faces().size());
    List<scalarList> posNOP(mesh.faces().size());

    wordList tableNames(thermo.composition().species());
    tableNames.append("he");
    Foam::combustionModels::tableSolver solver(Foam::combustionModels::tableSolver(mesh, tableNames));

    PtrList<volScalarField>& Y(thermo.composition().Y());
    volScalarField& he(thermo.he());

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << nl << "Time = " << runTime.timeName() << nl << endl;

        volScalarField Z
        (
            IOobject
            (
                "Z",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField chi
        (
            IOobject
            (
                "chi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        volScalarField varZ
        (
            IOobject
            (
                "varZ",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volScalarField Zeta
        (
            IOobject
            (
                "Zeta",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            sqrt(varZ/max(Z*(1-Z), SMALL))
       );
        
       volScalarField PVmajor
    (
        IOobject
        (
            "PVmajor",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );  //add by Jiang ÌíŒÓœø³Ì±äÁ¿
    volScalarField PVNO
    (
        IOobject
        (
            "PVNO",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );  //add by Jiang ÌíŒÓœø³Ì±äÁ¿

       // Interpolate for internal Field
       forAll(Y, i)
       {
    	  scalarField& YCells = Y[i].internalField();

          forAll(Z, cellI)
          {
        	 if (i == 0)
        	 {
                 if (useMixtureFractionVariance_) x[0] = min(Zeta[cellI], 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                 x[1] = Z[cellI];                                            //对于混合分数，直接将其赋值给起始边界网格
	   x[2] = PVmajor[cellI];                                            //对于进程变量，直接将其赋值给起始边界网格 	 
                 ubMajorIF[cellI] = solver.upperBounds(x);                                            //插值索引的上位点
                 posMajorIF[cellI] = solver.position(ubMajorIF[cellI], x);                                            //插值点的上位点及其值
	   
                if (useMixtureFractionVariance_) y[0] = min(Zeta[cellI], 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                 y[1] = Z[cellI];                                            //对于混合分数，直接将其赋值给起始边界网格
                 y[2] = PVNO[cellI];                                            //对于进程变量，直接将其赋值给起始边界网格      
                 ubNOIF[cellI] = solver.upperBounds(y);                                            //插值索引的上位点
                 posNOIF[cellI] = solver.position(ubNOIF[cellI], y);                                            //插值点的上位点及其值
        	 }
        	 // Info << tableNames[35]<<endl;
        	 if(i != 35)
        	 YCells[cellI] = solver.interpolate(ubMajorIF[cellI], posMajorIF[cellI], i);
              else
              YCells[cellI] = solver.interpolate(ubNOIF[cellI], posNOIF[cellI], i);           
          }
       }    

       forAll(he.boundaryField(), patchi)   
       {
          const fvPatchScalarField& pZ = Z.boundaryField()[patchi];
          const fvPatchScalarField& pPVmajor = PVmajor.boundaryField()[patchi];       //赋值进程变量的边界条件
          const fvPatchScalarField& pPVNO = PVNO.boundaryField()[patchi];       //赋值进程变量的边界条件	  
          fvPatchScalarField& pHe = he.boundaryField()[patchi];

          forAll(Y, i)
          {
        	  fvPatchScalarField& pY = Y[i].boundaryField()[patchi];

              forAll(pY , facei)//number of grids
              {
             	 if (i == 0)
             	 {
	       if (useMixtureFractionVariance_) x[0] = min(Zeta[facei], 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                     x[1] = pZ[facei];                                            //对于混合分数，直接将其赋值给起始边界网格
	       x[2] = pPVmajor[facei];                                            //对于进程变量，直接将其赋值给起始边界网格 

                    if (useMixtureFractionVariance_) y[0] = min(Zeta[facei], 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                    y[1] = pZ[facei];                                            //对于混合分数，直接将其赋值给起始边界网格
                    y[2] = pPVNO[facei];                                            //对于进程变量，直接将其赋值给起始边界网格 

                     ubMajorP[facei] = solver.upperBounds(x);
                     posMajorP[facei] = solver.position(ubMajorP[facei], x);      

                     ubNOP[facei] = solver.upperBounds(y);
                     posNOP[facei] = solver.position(ubNOP[facei], y);

                     pHe[facei] = solver.interpolate(ubMajorP[facei], posMajorP[facei], (solver.sizeTableNames() - 2));
             	 }
                  //Info << tableNames[35]<<endl;
                if(i != 35)
            	 pY[facei] = solver.interpolate(ubMajorP[facei], posMajorP[facei], i);
                else
              pY[facei] = solver.interpolate(ubNOP[facei], posNOP[facei], i);
             }
          }
       }   
        // Calculate thermodynamic Properties
        if (selectedFields.empty())
        {
        	forAll(Y, i)
            {
               Info << "Writing field " << thermo.composition().Y()[i].name() << endl;
        	   thermo.composition().Y()[i].write();
            }
            Info << "Writing field Zeta" << endl;
        	Zeta.write();
        }
        else
        {
        	forAll(Y, i)
            {
        	   if (selectedFields[thermo.composition().Y()[i].name()])
        	   {
                   Info << "Writing field " << thermo.composition().Y()[i].name() << endl;
        		   thermo.composition().Y()[i].write();
        	   }
            }
            Info << "Writing field Zeta" << endl;
        	Zeta.write();
        }

    }

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info<< "End\n" << endl;

	return 0;
}
