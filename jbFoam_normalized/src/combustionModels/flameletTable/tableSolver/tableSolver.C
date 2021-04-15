/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft

\*---------------------------------------------------------------------------*/


#include "tableSolver.H"

namespace Foam
{
namespace combustionModels
{

defineTypeNameAndDebug(tableSolver, 0);
defineRunTimeSelectionTable(tableSolver, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tableSolver::tableSolver(const fvMesh& mesh, const wordList& tableNames)
:
   IOdictionary
   (
      IOobject
      (
         "tableProperties",
         mesh.time().constant(),
         mesh,
         IOobject::MUST_READ_IF_MODIFIED,
         IOobject::NO_WRITE
      )
   ),
   tableNames_(tableNames),          //定义表的变量名数目
   paramNames_(3),          //定义插值索引的维度
   tables_(tableNames_.size() + 1),
   chiExt_(this->lookupOrDefault("chiExt", 0.0))
{
    forAll(tableNames_, i)
    {
        tableNames_[i] = tableNames_[i] + "_table";
        tables_.set(i, flameletTable::New(mesh, tableNames_[i]));
    }

    //paramNames_[0] = "chi_param";
    //paramNames_[1] = "Zeta_param";
    //paramNames_[2] = "Z_param";
    paramNames_[0] = "Zeta_param";
    paramNames_[1] = "Z_param";
    paramNames_[2] = "PV_param";   //第一维到第三维分别为混合分数的方差，混合分数，进程变量，para_为三线插值的基准坐标，基准空间

    forAll(paramNames_, i)
    {
    	params_.append(scalarList(this->lookup(paramNames_[i])));
    	if (paramNames_[i] == "chi_param" && chiExt_ != 0.0) params_[i].append(chiExt_);
    }                                  //给params_赋值
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tableSolver::~tableSolver()
{}

// * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * //

List<int> tableSolver::upperBounds(const scalarList& x) const
{
    // Upper, lower bound for table Interpolation and temp Value for bisection method
	List<int> ub(paramNames_.size(), 0);             //ub赋值为三维的列表
    int tlb, tub;
	int newVal = 0;

    // Determine upper bounds and interpolation weights for table interpolation
    // Bisection Method
    for (register label j=0; j<paramNames_.size(); j++)
    {
	   tub = params_[j].size() - 1;                          //将自定义的某一维列表最大值赋值给tub
   	   tlb = 0;
   	   while (tub-1 != tlb)
       {
          newVal = (tub + tlb)/2;
          if (x[j] < params_[j][newVal])
     	      tub = newVal;
          else
      	      tlb = newVal;
       }
   	   ub[j] = tub;
    }             //找到三个维度的上值点，赋值给ub。
   	// Return upper bounds
    return ub;
}

scalarList tableSolver::position(const List<int>& ub, const scalarList& x) const
{
    scalarList pos(paramNames_.size(), 0.0);

    for (register label j=0; j<paramNames_.size(); j++)
    {
        pos[j] = (x[j] - params_[j][ub[j]-1]) / (params_[j][ub[j]] - params_[j][ub[j]-1]);
    }

    return pos;
}             //三线差值的比例得出来

scalar tableSolver::interpolate(const List<int>& ub , const scalarList& pos, const label& i) const
{
    return tables_[i].interpolate(ub, pos);
}                    //对第i种变量进行插值

scalar tableSolver::PVcorrect(const List<int>& ub, const scalarList& pos, const label& i) const
{
    return tables_[i].PVcorrect(ub,pos);
}            

scalar tableSolver::maxChi()
{
        forAll(paramNames_, i)
        if (paramNames_[i] == "chi_param") return *params_[i].rbegin();
        return 0.0;
}

int tableSolver::sizeTableNames() const
{
        return tableNames_.size();
}

} // End combustionModels namespace
} // End Foam namespace
