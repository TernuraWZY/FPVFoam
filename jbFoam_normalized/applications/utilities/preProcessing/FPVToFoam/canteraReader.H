/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::canteraReader

Description
    Foam::canteraReader

SourceFiles
    canteraReader.C
    canteraLexer.C

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Gabriele Frank <gabriele.frank@unibw.de> Universität der Bundeswehr München

\*---------------------------------------------------------------------------*/

#ifndef canteraReader_H
#define canteraReader_H

#include "fileName.H"
#include "speciesTable.H"
#include "IOdictionary.H"
#include "List.H"
#include "HashTable.H"
#include "scalar.H"
#include "rhoReactionThermo.H"    
#include <FlexLexer.h>
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class chemkinReader Declaration
\*---------------------------------------------------------------------------*/

class canteraReader
:
    public yyFlexLexer
{

public:

private:

        static int yyBufSize;

    // Private Member Functions

        //- Flex lexer to read the cantera tables
        int lex();

        // Read the CHEMKIN files
        void read(const fileName& canteraFileName);

        void interpolateData();
		
	void interpolatePVMajorData();      //add by Jiang，用于将Chi维度转换为PV维度，并插值PV_param_的分布空间
	
	void interpolatePVNOData();		//added by wzy,读取NO的表

        void interpolateUnintegratedData();

        void calculateEnthalpy();

        void calculateThermo();

        void calculateZ();

        void betaPDFIntegration(const label& numChi, const scalar& varZ);

        //- Disallow default bitwise copy construct
        canteraReader(const canteraReader&);

        //- Disallow default bitwise assignment
        void operator=(const canteraReader&);

        //Private Data
        int col_iter;

        int num_lines, num_chars, num_columns;

        bool noSecondLine;

        basicMultiComponentMixture& composition;
        rhoReactionThermo& thermo;        

        PtrList<volScalarField>& Y_;
        //volScalarField& FI_;
        volScalarField& he_;             
        volScalarField& p_;

        hashedWordList  tableNames_;
        hashedWordList  tablesToBeRead_;

        List<scalar> 	chi_param_;     //保留这个表是为了读入用标量耗散率命名的文件
        List<scalar> 	Zeta_param_;
        List<scalar> 	Z_param_;
	List<scalar> 	PV_Major_param_;    //add by Jiang，进程变量的表。
	List<scalar> 	PV_NO_param_;    //added by wzy，进程变量NO的表。
	

        word 			mixtureFractionDefinition_;

        //Lists to manage the columns
        List<scalar> columns_;
        hashedWordList tableSorted_;

        //data for the interpolation routine
        List<scalar> coordinates_;
        int coord_col;
        scalar size_;

        List<List<List<List<scalar> > > > sampledMajorData_;
        List<List<List<List<scalar> > > > sampledNOData_;	//added by wzy,存储基于NO的使用表
	List<List<List<List<scalar> > > > sampledPVMajorData_;     //add by Jiang，用存储FPV使用表，及三维为混合分数方差，混合分数，及标量耗散率的表。
	List<List<List<List<scalar> > > > sampledPVNOData_;     //added by wzy，用存储NO使用表，及三维为混合分数方差，混合分数，及标量耗散率的表。

        List<scalar> enthalpyCantera_;
        int enthalpy_col;

        List<scalar> enthalpyOpenFoam_;
        List<scalar> enthalpyCMix_;
        List<List<scalar> > singleData_;
        List<List<scalar> > integratedData_;
        List<List<scalar> > unintegratedData_;

public:
        // Constructors

        //construct from filename
        canteraReader(const IOdictionary& canteraDict, rhoReactionThermo& thermo, basicMultiComponentMixture& composition);    

        //- Destructor
        virtual ~canteraReader();

        // Member functions
        int numberOfLines(){return num_lines;}

        int numberOfColumns(){return num_columns;}

        hashedWordList getNames();

        void write(const int& i,IOdictionary& dictionary,OFstream& output);
        
        void writeNO(const int& i,IOdictionary& dictionary,OFstream& output);
        
        void forewrite(const int& i,IOdictionary& dictionary,OFstream& output);
        
        void forewriteNO(const int& i,IOdictionary& dictionary,OFstream& output);
                
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
