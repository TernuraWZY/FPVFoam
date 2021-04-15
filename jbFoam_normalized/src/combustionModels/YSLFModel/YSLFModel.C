/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "YSLFModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType>
YSLFModel<CombThermoType>::YSLFModel
(
    const word& modelType, const fvMesh& mesh
)
:
    CombThermoType(modelType, mesh),
    solver_(tableSolver(mesh, tables())),
    Y_(this->thermo().composition().Y()),
    he_(this->thermo().he()), 
    T_(this->thermo().T()),	
                
    OMGmajor_(this->thermo().OMGmajor()),          //实例thermo的成员进程变量源项
    OMGNO_(this->thermo().OMGNO()),          //实例thermo的NO成员进程变量源项

    NO_(this->thermo().NO_Field()),
    NOTransport_(this->thermo().NOTransport()),
    OMGNOTransport_(this->thermo().OMGNOTransport()),
    OMGNOTransportPos_(this->thermo().OMGNOTransportPos()),
    OMGNOTransportNeg_(this->thermo().OMGNOTransportNeg()),
    PVmajor_(this->thermo().PVmajor()),          //实例thermo的成员进程变量
    PVNO_(this->thermo().PVNO()),          //实例thermo的NO成员进程变量    
    
    Z_(this->thermo().Z()),
    varZ_(this->thermo().varZ()),
    Chi_(this->thermo().Chi()),
        
    
    //FI_(this->thermo().FI()),
    ubmajorIF_(mesh.cells().size()),         //uppoint for mesh
    ubmajorP_(),
    ubNOIF_(mesh.cells().size()),         //uppoint for mesh
    ubNOP_(),                           //uppoint for face of mesh
                               //uppoint for face of mesh
    posmajorIF_(mesh.cells().size()),       //ratio for mesh
    posmajorP_(),                           //ratio for face of mesh     
    posNOIF_(mesh.cells().size()),       //ratio for mesh
    posNOP_(),                           //ratio for face of mesh
    
    useScalarDissipation_(this->coeffs().lookup("useScalarDissipation")),
    useMixtureFractionVariance_(this->coeffs().lookup("useMixtureFractionVariance"))
{
	const polyBoundaryMesh& patches = mesh.boundaryMesh();
	int patchSize = 0;
    forAll(patches, patchI)
    {
    	const polyPatch& pp = patches[patchI];
    	if (pp.size() > patchSize) patchSize = pp.size();
    }
	//Info<< "OMG is" <<OMGmajor_ <<endl;
    ubmajorP_.setSize(patchSize);
    posmajorP_.setSize(patchSize);
    ubNOP_.setSize(patchSize);
    posNOP_.setSize(patchSize);
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType>
YSLFModel<CombThermoType>::~YSLFModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CombThermoType>
hashedWordList YSLFModel<CombThermoType>::tables()
{
	hashedWordList tableNames = this->thermo().composition().species();
               tableNames.append("OMGNOTransport");
               tableNames.append("OMGNOTransportPos");
               tableNames.append("OMGNOTransportNeg");
	tableNames.append("PVmajor");
	tableNames.append("OMGmajor");
	tableNames.append("PVNO");
	tableNames.append("OMGNO");
	tableNames.append("he");
	tableNames.append("T");
	// for(int i = 0; i < tableNames.size(); i++)
	// Info<<tableNames[i]<<endl;
	return tableNames;
}

template<class CombThermoType>
void YSLFModel<CombThermoType>::correct()
{
    // limit the scalar dissipation rate to avoid instabilities at extinction
   // scalar chiLimiter = solver_.maxChi();
	Info<<"YSLF correct"<<endl;
    const scalarField& ZCells = Z_.internalField();
    const scalarField& varZCells = varZ_.internalField();

    const scalarField& PVmajorCells = PVmajor_.internalField();    //增加进程变量的赋值
    const scalarField& PVNOCells = PVNO_.internalField();    //增加进程变量的赋值

    scalarField& heCells = he_.internalField();
    scalarField& TCells = T_.internalField();
    scalarField& OMGmajorCells = OMGmajor_.internalField();       //增加进程变量源项的赋值
    scalarField& OMGNOCells = OMGNO_.internalField();       //增加进程变量源项的赋值
    scalarField& NOCells = NO_.internalField();
    scalarField& OMGNOTransportCells = OMGNOTransport_.internalField();
    scalarField& OMGNOTransportNegCells = OMGNOTransportNeg_.internalField();
    scalarField& OMGNOTransportPosCells = OMGNOTransportPos_.internalField();
    scalarField& NOTransportCells = NOTransport_.internalField();

    //- Update the species and enthalpy field
    if(this->active())
    {
       scalarList x(3, 0.0);//for major
       scalarList y(3, 0.0);//for NO
       double Zeta;

       // Interpolate for internal Field
       forAll(Y_, i)
       {
    	  scalarField& YCells = Y_[i].internalField();

          forAll(ZCells, cellI)
          {
        	 if (i == 0)
        	 {
        		 Zeta = sqrt(varZCells[cellI]/max(ZCells[cellI]*(1 - ZCells[cellI]), SMALL));     //用于限定varZ的范围
                 if (useMixtureFractionVariance_) x[0] = min(Zeta, 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                 x[1] = ZCells[cellI];                                          //对于混合分数，直接将其赋值给起始边界网格
		 x[2] = PVmajorCells[cellI];                                            //对于进程变量，直接将其赋值给起始边界网格 	
		 	 
		 y[0] = min(Zeta, 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                 y[1] = ZCells[cellI];                                            //对于混合分数，直接将其赋值给起始边界网格
		 y[2] = PVNOCells[cellI];   
		 ubNOIF_[cellI] = solver_.upperBounds(y);                                            //插值索引的上位点             
                 posNOIF_[cellI] = solver_.position(ubNOIF_[cellI], y);                                            //插值点的上位点及其值	 
                 ubmajorIF_[cellI] = solver_.upperBounds(x);                                            //插值索引的上位点
                 posmajorIF_[cellI] = solver_.position(ubmajorIF_[cellI], x);                                            //插值点的上位点及其值
                 
		 OMGmajorCells[cellI] = solver_.interpolate(ubmajorIF_[cellI], posmajorIF_[cellI], (solver_.sizeTableNames() - 5));    
		 OMGNOCells[cellI] = solver_.interpolate(ubNOIF_[cellI], posNOIF_[cellI], (solver_.sizeTableNames() - 3));
                             OMGNOTransportPosCells[cellI] = solver_.interpolate(ubNOIF_[cellI], posNOIF_[cellI], (solver_.sizeTableNames() - 8));
                             OMGNOTransportNegCells[cellI] = solver_.interpolate(ubNOIF_[cellI], posNOIF_[cellI], (solver_.sizeTableNames() - 7));
		 
                 heCells[cellI] = solver_.interpolate(ubmajorIF_[cellI], posmajorIF_[cellI], (solver_.sizeTableNames() - 2));       
                 TCells[cellI] = solver_.interpolate(ubmajorIF_[cellI], posmajorIF_[cellI], (solver_.sizeTableNames() - 1));                                   	 
        	 }
               if(i != 35 && i != 36 && i != 37)       	        	 
        	 YCells[cellI] = solver_.interpolate(ubmajorIF_[cellI], posmajorIF_[cellI], i); 
             else
                YCells[cellI] = solver_.interpolate(ubNOIF_[cellI], posNOIF_[cellI], i);            
            if(i == 35)
            {
                NOCells[cellI] = YCells[cellI];
             //   Info << "SMALL"<<endl;
                OMGNOTransportCells[cellI] = 3*OMGNOTransportNegCells[cellI]* ((NOTransportCells[cellI])/(NOCells[cellI] ))+OMGNOTransportPosCells[cellI];
              //  Info<<"SMALL successs"<<endl;
            }
          }
       }    

       forAll(he_.boundaryField(), patchi)   
       {
          const fvPatchScalarField& pvarZ = varZ_.boundaryField()[patchi];
          const fvPatchScalarField& pZ = Z_.boundaryField()[patchi];
          const fvPatchScalarField& pPVmajor = PVmajor_.boundaryField()[patchi];       //赋值进程变量的边界条件
          const fvPatchScalarField& pPVNO = PVNO_.boundaryField()[patchi];
		  
          fvPatchScalarField& pHe = he_.boundaryField()[patchi];
          fvPatchScalarField& pT = T_.boundaryField()[patchi];
          fvPatchScalarField& pOMGmajor = OMGmajor_.boundaryField()[patchi];       //赋值进程变量源项的边界条件
          fvPatchScalarField& pOMGNO = OMGNO_.boundaryField()[patchi];
          fvPatchScalarField& pNO = NO_.boundaryField()[patchi];
    fvPatchScalarField& pOMGNOTransport = OMGNOTransport_.boundaryField()[patchi];
    fvPatchScalarField& pOMGNOTransportNeg = OMGNOTransportNeg_.boundaryField()[patchi];
    fvPatchScalarField& pOMGNOTransportPos = OMGNOTransportPos_.boundaryField()[patchi];
    fvPatchScalarField& pNOTransport = NOTransport_.boundaryField()[patchi];
          forAll(Y_, i)
          {
        	  fvPatchScalarField& pY = Y_[i].boundaryField()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                     Zeta = sqrt(pvarZ[facei]/max(pZ[facei]*(1 - pZ[facei]), SMALL));
		    if (useMixtureFractionVariance_) x[0] = min(Zeta, 0.99);         //限定varZ的范围不超过Z(Z-1)，起始边界网格的方差不超过在Z（Z-1）
                     x[1] = pZ[facei];                                            //对于混合分数，直接将其赋值给起始边界网格
		     x[2] = pPVmajor[facei];                                            //对于进程变量，直接将其赋值给起始边界网格 
		     if (useMixtureFractionVariance_) y[0] = min(Zeta, 0.99); 
		     y[1] = pZ[facei];                                            //对于混合分数，直接将其赋值给起始边界网格NO
		     y[2] = pPVNO[facei];                                            //对于进程变量，直接将其赋值给起始边界网格NO

                     ubmajorP_[facei] = solver_.upperBounds(x);
                     posmajorP_[facei] = solver_.position(ubmajorP_[facei], x); 
                     
                     ubNOP_[facei] = solver_.upperBounds(y);
                     posNOP_[facei] = solver_.position(ubNOP_[facei], y);
                   
                     pOMGmajor[facei] = solver_.interpolate(ubmajorP_[facei], posmajorP_[facei], (solver_.sizeTableNames() - 5));                                         
                     pHe[facei] = solver_.interpolate(ubmajorP_[facei], posmajorP_[facei], (solver_.sizeTableNames() - 2));

		     pOMGNO[facei] = solver_.interpolate(ubNOP_[facei], posNOP_[facei], (solver_.sizeTableNames() - 3));

                             pOMGNOTransportPos[facei] = solver_.interpolate(ubNOP_[facei], posNOP_[facei], (solver_.sizeTableNames() - 8));
                             pOMGNOTransportNeg[facei] = solver_.interpolate(ubNOP_[facei], posNOP_[facei], (solver_.sizeTableNames() - 7));

		     pT[facei] = solver_.interpolate(ubmajorP_[facei], posmajorP_[facei], (solver_.sizeTableNames() - 1));
             	 }
                 if(i != 35 && i != 36 && i != 37)
            	 pY[facei] = solver_.interpolate(ubmajorP_[facei], posmajorP_[facei], i); 
                 else
                pY[facei] = solver_.interpolate(ubNOP_[facei], posNOP_[facei], i);   
                if(i == 35)
                {
                    pNO[facei] = pY[facei];
                    pOMGNOTransport[facei] = 3*pOMGNOTransportNeg[facei]* ((pNOTransport[facei] ) /(pNO[facei] )) + pOMGNOTransportPos[facei];
                }
             }
          }
       }
       this->thermo().correct();
    }
}

template<class CombThermoType>
Switch YSLFModel<CombThermoType>::correctDensity()
{
	return true;
}

template<class CombThermoType>
Foam::tmp<Foam::fvScalarMatrix>
YSLFModel<CombThermoType>::R
(
    volScalarField& Y              
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    return tSu;
}

template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
YSLFModel< CombThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tSh;
}

template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
YSLFModel< CombThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tdQ;
}

template<class CombThermoType>
bool YSLFModel<CombThermoType>::read()
{
    if (CombThermoType::read())
    {
        this->coeffs().lookup("useScalarDissipation") >> useScalarDissipation_;
        this->coeffs().lookup("useMixtureFractionVariance") >> useMixtureFractionVariance_;
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
