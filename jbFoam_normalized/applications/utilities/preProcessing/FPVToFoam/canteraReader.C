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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Gabriele Frank <gabriele.frank@unibw.de> Universität der Bundeswehr München

\*---------------------------------------------------------------------------*/
#include "canteraReader.H"
#include "IFstream.H"
#include "fileName.H"
#include "OFstream.H"
#include "IOList.H"
#include "Gamma.h"
#include <limits>
#include <sstream>
#include <fstream>

void Foam::canteraReader::read(const fileName& canteraFileName)
{
    //clear the data containers
    singleData_.clear();
    singleData_.setSize(tablesToBeRead_.size());
    tableSorted_.clear();
    coordinates_.clear();
    enthalpyCantera_.clear();

    std::ifstream myinfile(canteraFileName.c_str());
    yy_buffer_state* bufferPtr(yy_create_buffer(&myinfile, yyBufSize));
    yy_switch_to_buffer(bufferPtr);
    col_iter=0;
    num_columns=0;
    num_lines=0;
    noSecondLine=1;
    while (lex()!=0)
    {}
    yy_delete_buffer(bufferPtr);
}

void Foam::canteraReader::betaPDFIntegration(const label& numChi, const scalar& Zeta)
{
if (Zeta != 0)
{
   List<scalar> Z_(integratedData_[tableNames_["Z"]]);
   List<scalar> varZ_(integratedData_[tableNames_["Z"]].size(), 0.0);
   List<scalar> pdfAlpha(integratedData_[tableNames_["Z"]].size(), 0.0);
   List<scalar> pdfBeta(integratedData_[tableNames_["Z"]].size(), 0.0);
   List<scalar> PDF(integratedData_[tableNames_["Z"]].size(), 0.0);

   for (int i=0;i<integratedData_[tableNames_["Z"]].size();i++)
   {
      varZ_[i] = sqr(Zeta) * (Z_[i]*(1.0 - Z_[i]));

      if (varZ_[i] > 1e-7)
      {
          pdfAlpha[i] = Z_[i] * ((Z_[i]*(1.0-Z_[i]))/varZ_[i]-1);
          pdfBeta[i] = (1.0 - Z_[i]) * ((Z_[i]*(1.0-Z_[i]))/varZ_[i]-1);

          // Limit alpha and beta but keep their ratio
          if (pdfAlpha[i] > 500)
          {
             pdfBeta[i] = pdfBeta[i]/pdfAlpha[i] * 500;
             pdfAlpha[i] = 500;
          }
          if (pdfBeta[i] > 500)
          {
             pdfAlpha[i] = pdfAlpha[i]/pdfBeta[i] * 500;
             pdfBeta[i] = 500;
          }

          int    gridPoints = 250;
          List<long double> hZ_(gridPoints, 0.0);
          List<long double> helpZ_(gridPoints, 0.0);
          List<long double> delta_(gridPoints, 0.0);

          if ((pdfAlpha[i] > 1) && (pdfBeta[i] > 1))
          {
             // Allocation of Z for PDF integration
             scalar Zmax = 0;
             int    n1 = 0;
             int    n2 = 0;
             PDF.clear();
             PDF.resize(Z_.size(), 0.0);
             
             for (int j=0;j<Z_.size();j++)
             {
                PDF[j] = std::pow(Z_[j],(pdfAlpha[i]-1.0)) * std::pow((1.0 - Z_[j]),(pdfBeta[i]-1.0)) * Gamma(pdfAlpha[i] + pdfBeta[i]) / (min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
                if (PDF[j] > PDF[max(j-1, 0)]) Zmax = Z_[j];
             }

             if(pdfAlpha[i]/pdfBeta[i] <= 0.5)
             {
                n1 = 0.2*gridPoints;
                n2 = 0.8*gridPoints+1;
             }
             else if (pdfAlpha[i]/pdfBeta[i] >= 2)
             {
                n1 = 0.8*gridPoints;
                n2 = 0.2*gridPoints+1;
             }
             else
             {
                n1 = 0.5*gridPoints;
                n2 = 0.5*gridPoints+1;
             }

             //  Allocate Z for 0 < Z < Zmax
             scalar ex1 = 0.9;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             for (int j=0; j<n1; j++)
             {
                delta_[j] = pow(ex1,j) * delta_[0];
                hZ_[j+1] = hZ_[j] + delta_[j];
             }
             for (int j=1; j<n1; j++)
             {
                hZ_[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 1.1;
             delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
               delta_[j] = pow(ex2,j) * delta_[0];
               helpZ_[j+1] = helpZ_[j] + delta_[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ_[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hZ_[j] /= hZ_[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hZ_.size(), 0.0);

             for (int j=0;j<hZ_.size();j++)
             {
                PDF[j] = (std::pow(hZ_[j],(pdfAlpha[i]-1e0))) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1e0)) * Gamma(pdfAlpha[i] + pdfBeta[i])/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
             }
          }

          else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] > 1))
          {
             // PDF Singularity at Z = 0
             // Allocation of Z for PDF integration
             int    n1 = 0;
             int    n2 = 0;

             if (pdfAlpha[i]/pdfBeta[i] > 0.5)
             {
                scalar Zmax = 0.5;
                scalar ex1 = 1.1;
                n1 = 0.7*gridPoints;
                delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

                // Allocate Z for 0 < Z < Zmax
                for (int j=0; j<n1; j++)
                {
                   delta_[j] = pow(ex1,j) * delta_[0];
                   hZ_[j+1] = hZ_[j] + delta_[j];
                }
                for (int j=1; j<n1; j++)
                {
                   hZ_[j] *= Zmax;
                }

                // Allocate Z for Zmax < Z < 1
                scalar ex2 = 1.1;
                n2 = 0.3*gridPoints+1;
                delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

                for (int j=0; j<n2-1; j++)
                {
                   delta_[j] = pow(ex2,j) * delta_[0];
                   helpZ_[j+1] = helpZ_[j] + delta_[j];
                }
                for (int j=0;j<n2;j++)
                {
                   helpZ_[j] *= (1.0 - Zmax);
                }
                for (int j=0;j<n2-1;j++)
                {
                   hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
                }
             }

             else
	         {
                scalar ex2 = 1.05;
                delta_[0] = (1.0 - ex2)/(1.0 - pow(ex2,(gridPoints-1)));
                for (int j=0; j<gridPoints-1; j++)
                {
                   delta_[j] = pow(ex2,j) * delta_[0];
                   hZ_[j+1] = hZ_[j] + delta_[j];
                }
             }

             // Scaling
             for (int j=0;j<gridPoints;j++)
             {
                hZ_[j] /= hZ_[gridPoints-1];
             }

             // Calculate BetaPDF
             PDF.clear();
             PDF.resize(hZ_.size(), 0.0);
             for (int j=1;j<hZ_.size();j++)
             {
                PDF[j] = std::pow(hZ_[j],(pdfAlpha[i]-1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
             }
             PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];

          }    

          else if ((pdfAlpha[i] > 1) && (pdfBeta[i] <= 1))
          {
          // PDF Singularity at Z = 1
          // Allocation of Z for PDF integration

          int    n1 = 0;
          int    n2 = 0;

          if (pdfAlpha[i]/pdfBeta[i] < 2)
          {
             scalar Zmax = 0.5;
             scalar ex1 = 1.1;
             n1 = 0.3*gridPoints;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

             // Allocate Z for 0 < Z < Zmax
             for (int j=0; j<n1; j++)
             {
                delta_[j] = pow(ex1,j) * delta_[0];
                hZ_[j+1] = hZ_[j] + delta_[j];
             }
             for (int j=1; j<n1; j++)
             {
                hZ_[j] *= Zmax;
             }

             // Allocate Z for Zmax < Z < 1
             scalar ex2 = 0.9;
             n2 = 0.7*gridPoints+1;
             delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

             for (int j=0; j<n2-1; j++)
             {
                delta_[j] = pow(ex2,j) * delta_[0];
                helpZ_[j+1] = helpZ_[j] + delta_[j];
             }
             for (int j=0;j<n2;j++)
             {
                helpZ_[j] *= (1.0 - Zmax);
             }
             for (int j=0;j<n2-1;j++)
             {
                hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
             }
          }

          else
          {
             scalar ex1 = 0.95;
             delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(gridPoints-1)));

             for (int j=0; j<gridPoints-1; j++)
             {
                 delta_[j] = pow(ex1,j) * delta_[0];
                 hZ_[j+1] = hZ_[j] + delta_[j];
             }
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hZ_[j] /= hZ_[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hZ_.size(), 0.0);

          for (int j=0;j<hZ_.size()-1;j++)
          {
              PDF[j] = std::pow((hZ_[j]),(pdfAlpha[i]- 1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];

       }

       else if ((pdfAlpha[i] <= 1) && (pdfBeta[i] <= 1))
       {
          // PDF Singularity at Z = 1 and Z = 0
          // Allocation of Z for PDF integration
          int    n1 = 0;
          int    n2 = 0;
          scalar Zmax = 0.5;
          scalar ex1 = 1.1;
          n1 = 0.5*gridPoints;
          delta_[0] = (1.0 - ex1)/(1.0 - pow(ex1,(n1-1)));

          // Allocate Z for 0 < Z < Zmax
          for (int j=0; j<n1; j++)
          {
              delta_[j] = pow(ex1,j) * delta_[0];
              hZ_[j+1] = hZ_[j] + delta_[j];
          }
          for (int j=1; j<n1; j++)
          {
             hZ_[j] *= Zmax;
          }

          // Allocate Z for Zmax < Z < 1
          scalar ex2 = 0.9;
          n2 = 0.5*gridPoints+1;
          delta_[0] = (1.0-ex2)/(1.0-pow(ex2,(n2-1)));

          for (int j=0; j<n2-1; j++)
          {
              delta_[j] = pow(ex2,j) * delta_[0];
              helpZ_[j+1] = helpZ_[j] + delta_[j];
          }
          for (int j=0;j<n2;j++)
          {
              helpZ_[j] *= (1.0 - Zmax);
          }
          for (int j=0;j<n2-1;j++)
          {
             hZ_[n1+j] = hZ_[n1-1]+helpZ_[j];
          }

          // Scaling
          for (int j=0;j<gridPoints;j++)
          {
             hZ_[j] /= hZ_[gridPoints-1];
          }

          // Calculate BetaPDF
          PDF.clear();
          PDF.resize(hZ_.size(), 0.0);

          for (int j=1;j<hZ_.size()-1;j++)
          {
             PDF[j] = std::pow(hZ_[j],(pdfAlpha[i]- 1e0)) * std::pow((1e0 - hZ_[j]),(pdfBeta[i]-1.0)) * min(Gamma(pdfAlpha[i] + pdfBeta[i]),1e17)/(min(Gamma(pdfAlpha[i]),1e17)*min(Gamma(pdfBeta[i]),1e17));
          }
          PDF[gridPoints - 1] = 1.5 * PDF[gridPoints - 2] / pdfBeta[i];
          PDF[0] = 1.5 * PDF[1] / pdfAlpha[i];
       }

       // Calculate the area of the PDF for scaling
       scalar intPDF = 0;
       for (int j=1;j<hZ_.size();j++)
       {
    	  intPDF += (hZ_[j-1] - hZ_[j]) * (PDF[j-1] + PDF[j])/2;
       }

       // Interpolate singleData entries to the new mixture fraction space
       List<scalar> hY_(gridPoints, 0.0);
       scalar intY = 0;

       for (int j=0;j<integratedData_.size();j++)
       {
          hY_ = 0.0;
          hY_[0] = singleData_[j][0];
          hY_[hY_.size()-1] = singleData_[j][singleData_[j].size()-1];
          intY = 0;

          for (int k=1;k<hZ_.size()-1;k++)
          {
             int ubZ = 0;
             for (int l=0;l<Z_.size();l++)
             {
                ubZ = l;
                if (hZ_[k] < Z_[l])
                break;
             }
             int lbZ = ubZ -1;

             // Interpolation to hZ space
             hY_[k] = (singleData_[j][ubZ] - singleData_[j][lbZ])/max(Z_[ubZ] - Z_[lbZ], SMALL) * (hZ_[k] - Z_[lbZ]) + singleData_[j][lbZ];
             // PDF Integration using the trapezoidal rule
             intY += (hZ_[k-1] - hZ_[k]) * (hY_[k-1]*PDF[k-1] + hY_[k]*PDF[k])/(2.0 * intPDF);
          }

          // Special treatment for the boundaries
          intY += (hZ_[hZ_.size()-2] - hZ_[hZ_.size()-1]) * (hY_[hZ_.size()-2]*PDF[hZ_.size()-2] + hY_[hZ_.size()-1]*PDF[hZ_.size()-1])/(2.0 * intPDF);
          if (i != 0 && i != integratedData_[tableNames_["Z"]].size()-1 && tableNames_[j] != "Z")
             integratedData_[j][i] = intY;
       }
     }
   }
 }
}

void Foam::canteraReader::interpolateData()
{
List<scalar> Z_(integratedData_[tableNames_["Z"]]);
List<scalar> hY_(Z_param_.size(), 0.0);

for (int i=0;i<tableNames_.size();i++)
{
   hY_ = 0;
   for (int j=0;j<Z_param_.size();j++)
   {
      int ubZ = 0;
      for (int k=0;k<Z_.size();k++)
      {
         ubZ = k;
         if (Z_[k] > Z_param_[j])
            break;
      }
      int lbZ = ubZ - 1;
      hY_[j] = (integratedData_[i][ubZ] - integratedData_[i][lbZ])/max(Z_[ubZ] - Z_[lbZ], SMALL) * (Z_param_[j] - Z_[lbZ]) + integratedData_[i][lbZ];
    }
    hY_[0] = integratedData_[i][0];
    hY_[Z_param_.size()-1] = integratedData_[i][Z_.size()-1];
    integratedData_[i] = hY_;
}
}



void Foam::canteraReader::interpolatePVMajorData()      //add by jiang 换维函数，将标量耗散率维度换成进程变量维度。
{
  //List<scalar> PV_(sampledMajorData_[tableNames_["PV"]]);
  List<scalar> hY_(PV_Major_param_.size(), 0.0);
  Info << "init hY size is " << hY_.size() << endl;

   for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
   {
     for (int numZ=0; numZ<Z_param_.size();numZ++)
     {
        double temp=0.00;
        for(int j = 0;j < chi_param_.size();j++)
        {
        
            for(int k=chi_param_.size()-1;k>j;k--)
            {
            //Info << "swap process" << endl;
                if(sampledMajorData_[tableNames_["PVmajor"]][k][numZeta][numZ] < sampledMajorData_[tableNames_["PVmajor"]][k-1][numZeta][numZ])
                {
		            for (int i=0;i<tableNames_.size();i++)
                    {
                       temp =sampledMajorData_[i][k][numZeta][numZ];
                       sampledMajorData_[i][k][numZeta][numZ] = sampledMajorData_[i][k-1][numZeta][numZ];
                       sampledMajorData_[i][k-1][numZeta][numZ] = temp;
                    }
                }
	    }
        }
        //Info << "chi_param size is " << chi_param_.size() <<endl;
        if(numZ!=Z_param_.size()-1 && numZ!=0)
          for (int j=0;j<chi_param_.size();j++)
             {
		//Info << "normalize  " << j << endl;
		if (sampledMajorData_[tableNames_["PVmajor"]][chi_param_.size()-1][numZeta][numZ]<1e-7)
		{
		//Info << "the divisor is " << sampledMajorData_[tableNames_["PVmajor"]][chi_param_.size()-1][numZeta][numZ] <<endl;	
		}
		
                sampledMajorData_[tableNames_["PVmajor"]][j][numZeta][numZ] /=  sampledMajorData_[tableNames_["PVmajor"]][chi_param_.size()-1][numZeta][numZ];
               // Info << "normalize OMG count: " << j << endl;
                sampledMajorData_[tableNames_["OMGmajor"]][j][numZeta][numZ] /=  sampledMajorData_[tableNames_["PVmajor"]][chi_param_.size()-1][numZeta][numZ];
               //Info << "total count: " << j << endl; 
             }
      }
      Info << "normalize over" << endl;
   }  //用冒泡法将数组按chi从小到大变成按PV从小到大的顺序排序
   Info << "from big to small" << endl;
    for (int i=0;i<tableNames_.size();i++)
    {
        for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
        {
            //for (int numZ=1; numZ<Z_param_.size()-1;numZ++)
            for (int numZ=0; numZ<Z_param_.size();numZ++)
            {      
            
               for (int numPV=1; numPV<PV_Major_param_.size()-1;numPV++)
                {
               // Info<< "numZ is "<<numZ<<" and PV_param_  "<<numPV<< " is " <<PV_Major_param_[numPV]<<endl;
                 if (numZ==0)
                   hY_[numPV]=sampledMajorData_[i][0][numZeta][0];
                 else if (numZ==Z_param_.size()-1)
                   hY_[numPV]=sampledMajorData_[i][0][numZeta][Z_param_.size()-1];	     
                 else
                 {           
	          int ubPV = 0;
		  for (int k=0; k<chi_param_.size();k++)
	           {
                    ubPV=k;
                   // Info<< "sampledMajorData_ is "<<sampledMajorData_[tableNames_["PVmajor"]][ubPV][numZeta][numZ]
                   // <<" PV_Major_param_[numPV] is "<<PV_Major_param_[numPV]<<endl;
		    if (sampledMajorData_[tableNames_["PVmajor"]][ubPV][numZeta][numZ] > PV_Major_param_[numPV])
                       break; 
	           }         
                  int lbPV = ubPV - 1;
//                  if(lbPV < 0)
//                  break;
                //  Info<<"lbPV is "<<lbPV<<" and ubPV is "<<ubPV <<endl;
                  hY_[numPV] = (sampledMajorData_[i][ubPV][numZeta][numZ] - sampledMajorData_[i][lbPV][numZeta][numZ])/max(sampledMajorData_[tableNames_["PVmajor"]][ubPV][numZeta][numZ] - sampledMajorData_[tableNames_["PVmajor"]][lbPV][numZeta][numZ], SMALL) * (PV_Major_param_[numPV] - sampledMajorData_[tableNames_["PVmajor"]][lbPV][numZeta][numZ]) + sampledMajorData_[i][lbPV][numZeta][numZ];
                  //Info<<"interpolate over" <<endl;
                  }				  
	        }
                hY_[0]=sampledMajorData_[i][0][numZeta][numZ];
                hY_[PV_Major_param_.size()-1]=sampledMajorData_[i][chi_param_.size()-1][numZeta][numZ];
		sampledPVMajorData_[i][numZeta][numZ] = hY_; 
            }     //将标量耗散率、混合分数方差、混合分数为索引的表转化为以混合分数方差、混合分数、进程变量为索引的表。
	    sampledPVMajorData_[i][numZeta][0][0]=sampledMajorData_[i][0][numZeta][0];
            sampledPVMajorData_[i][numZeta][Z_param_.size()-1][0]=sampledMajorData_[i][0][numZeta][Z_param_.size()-1];            
        }       		    
    }   
   // Info << "interpolatePVMajorData over" << endl;
}


void Foam::canteraReader::interpolatePVNOData()      //add by jiang 换维函数，将标量耗散率维度换成进程变量维度。
{
  //List<scalar> PV_(sampledMajorData_[tableNames_["PV"]]);
  List<scalar> hY_(PV_NO_param_.size(), 0.0);
  Info << "init hY size is " << hY_.size() << endl;

   for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
   {
     for (int numZ=0; numZ<Z_param_.size();numZ++)
     {
        double temp=0.00;
        for(int j = 0;j < chi_param_.size();j++)
        {
        
            for(int k=chi_param_.size()-1;k>j;k--)
            {
            //Info << "swap process" << endl;
                if(sampledNOData_[tableNames_["PVNO"]][k][numZeta][numZ] < sampledNOData_[tableNames_["PVNO"]][k-1][numZeta][numZ])
                {
		            for (int i=0;i<tableNames_.size();i++)
                    {
                       temp =sampledNOData_[i][k][numZeta][numZ];
                       sampledNOData_[i][k][numZeta][numZ] = sampledNOData_[i][k-1][numZeta][numZ];
                       sampledNOData_[i][k-1][numZeta][numZ] = temp;
                    }
                }
	    }
        }
        //Info << "chi_param size is " << chi_param_.size() <<endl;
        if(numZ!=Z_param_.size()-1 && numZ!=0)
          for (int j=0;j<chi_param_.size();j++)
             {
		//Info << "normalize  " << j << endl;
		if (sampledNOData_[tableNames_["PVNO"]][chi_param_.size()-1][numZeta][numZ]<1e-7)
		{
		//Info << "the divisor is " << sampledMajorData_[tableNames_["PVmajor"]][chi_param_.size()-1][numZeta][numZ] <<endl;	
		}
		
                sampledNOData_[tableNames_["PVNO"]][j][numZeta][numZ] /=  sampledNOData_[tableNames_["PVNO"]][chi_param_.size()-1][numZeta][numZ];
               // Info << "normalize OMG count: " << j << endl;
                sampledNOData_[tableNames_["OMGNO"]][j][numZeta][numZ] /=  sampledNOData_[tableNames_["PVNO"]][chi_param_.size()-1][numZeta][numZ];
               //Info << "total count: " << j << endl; 
             }
      }
      Info << "normalize over" << endl;
   }  //用冒泡法将数组按chi从小到大变成按PV从小到大的顺序排序
   Info << "from big to small" << endl;
    for (int i=0;i<tableNames_.size();i++)
    {
        for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
        {
            //for (int numZ=1; numZ<Z_param_.size()-1;numZ++)
            for (int numZ=0; numZ<Z_param_.size();numZ++)
            {      
            
               for (int numPV=1; numPV<PV_NO_param_.size()-1;numPV++)
                {
               // Info<< "numZ is "<<numZ<<" and PV_param_  "<<numPV<< " is " <<PV_Major_param_[numPV]<<endl;
                 if (numZ==0)
                   hY_[numPV]=sampledNOData_[i][0][numZeta][0];
                 else if (numZ==Z_param_.size()-1)
                   hY_[numPV]=sampledNOData_[i][0][numZeta][Z_param_.size()-1];	     
                 else
                 {           
	          int ubPV = 0;
		  for (int k=0; k<chi_param_.size();k++)
	           {
                    ubPV=k;
		    if (sampledNOData_[tableNames_["PVNO"]][ubPV][numZeta][numZ] > PV_NO_param_[numPV])
                       break; 
	           }         
                  int lbPV = ubPV - 1;
                  hY_[numPV] = (sampledNOData_[i][ubPV][numZeta][numZ] - sampledNOData_[i][lbPV][numZeta][numZ])/max(sampledNOData_[tableNames_["PVNO"]][ubPV][numZeta][numZ] - sampledNOData_[tableNames_["PVNO"]][lbPV][numZeta][numZ], SMALL) * (PV_NO_param_[numPV] - sampledNOData_[tableNames_["PVNO"]][lbPV][numZeta][numZ]) + sampledNOData_[i][lbPV][numZeta][numZ];
                  //Info<<"interpolate over" <<endl;
                  }				  
	        }
                hY_[0]=sampledNOData_[i][0][numZeta][numZ];
                hY_[PV_NO_param_.size()-1]=sampledNOData_[i][chi_param_.size()-1][numZeta][numZ];
		sampledPVNOData_[i][numZeta][numZ] = hY_; 
            }     //将标量耗散率、混合分数方差、混合分数为索引的表转化为以混合分数方差、混合分数、进程变量为索引的表。
	    sampledPVNOData_[i][numZeta][0][0]=sampledNOData_[i][0][numZeta][0];
            sampledPVNOData_[i][numZeta][Z_param_.size()-1][0]=sampledNOData_[i][0][numZeta][Z_param_.size()-1];            
        }       		    
    }   
   // Info << "interpolatePVMajorData over" << endl;
}

void Foam::canteraReader::calculateEnthalpy()
{
Info << "calculate sensible Enthalpy" << endl;

scalar pstd = 1.01325e5;
label labelT(tableNames_["T"]);
List<scalar> he(singleData_[labelT].size(), 0.0);

for (int i=0;i<singleData_[labelT].size();i++)
{
   for (int j=0; j<thermo.composition().species().size();j++)
   {
      label k = composition.species()[tableNames_[j]];
      he[i] += singleData_[j][i] * composition.Hs(k, pstd, singleData_[labelT][i]);
   }
}
singleData_.append(he);
}

void Foam::canteraReader::calculateZ()
{

if (mixtureFractionDefinition_ == "readFromTable")
{
   Info << "read mixture Fraction from Table" << endl;
}

else
{
   FatalErrorIn("Foam::canteraReader::calculateZ()")
   << "Unknown mixture fraction definition " << nl
   << "Valid mixture fraction definition are :" << nl
   << "readFromTable"  << nl
   << exit(FatalIOError);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::canteraReader::canteraReader(const IOdictionary& canteraDict, rhoReactionThermo& thermo, basicMultiComponentMixture& composition)  
:  composition(composition),
   thermo(thermo),
   Y_(thermo.composition().Y()),
   he_(thermo.he()),      
   p_(thermo.p()),
   tableNames_(thermo.composition().species()),
   tablesToBeRead_(thermo.composition().species()),
   chi_param_(canteraDict.lookup("chi_param")),
   Zeta_param_(canteraDict.lookup("Zeta_param")),
   Z_param_(canteraDict.lookup("Z_param")),
   PV_Major_param_(canteraDict.lookup("PV_param")),   //add by jiang，进程变量的索引分布
   PV_NO_param_(canteraDict.lookup("PV_param")),	//added by wzy,读取NO表
   mixtureFractionDefinition_(canteraDict.lookup("mixtureFractionDefinition")),
   columns_(tableNames_.size()+10)
{
   p_ = dimensionedScalar("p", dimPressure, canteraDict.lookup("operatingPressure"));

   tablesToBeRead_.append("T");
   tablesToBeRead_.append("Z");
   tablesToBeRead_.append("PVmajor");    //add by jiang, read PV
   tablesToBeRead_.append("OMGmajor");   //add by jiang, read OMG
   tablesToBeRead_.append("PVNO");    //add by wzy, read PVNO
   tablesToBeRead_.append("OMGNO");   //add by wzy, read OMGNO
  // tablesToBeRead_.append("NOTransport");    //add by wzy, 
   //tablesToBeRead_.append("OMGNOTransport");   //add by wzy, read OMGNO
   tablesToBeRead_.append("OMGNOTransportPos");   //add by wzy, read OMGNO
   tablesToBeRead_.append("OMGNOTransportNeg");   //add by wzy, read OMGNO

   tableNames_.append("T");
   tableNames_.append("Z");
   //tableNames_.append("he");
   tableNames_.append("PVmajor");       //add by jiang, add PV into table
   tableNames_.append("OMGmajor");       //add by jiang, add OMG into table
   tableNames_.append("PVNO");       //add by wzy, add PV into table
   tableNames_.append("OMGNO");       //add by wzy, add OMG into table
   //tableNames_.append("NOTransport");    //add by wzy, 
  // tableNames_.append("OMGNOTransport");   //add by wzy, read OMGNO
   tableNames_.append("OMGNOTransportPos");   //add by wzy, read OMGNO
   tableNames_.append("OMGNOTransportNeg");   //add by wzy, read OMGNO
   tableNames_.append("he");
   
   //读取maxNO，来进行NO非均一化选点 added by wzy  
   
   fileName canteraFileName(canteraDict.lookup("canteraFileName"));

   //sampledData_ correct size of the Lists
   sampledMajorData_.resize(tableNames_.size());

   for (int i=0; i<sampledMajorData_.size(); i++)
   {
      sampledMajorData_[i].resize(chi_param_.size());
      for (int j=0; j<sampledMajorData_[i].size();j++)
      {
         sampledMajorData_[i][j].resize(Zeta_param_.size());
         for (int k=0; k<sampledMajorData_[i][j].size();k++)
         {
            sampledMajorData_[i][j][k].resize(Z_param_.size());
         }
      }
   }
   //added by wzy 存NO数据的表
   sampledNOData_.resize(tableNames_.size());

   for (int i=0; i<sampledNOData_.size(); i++)
   {
      sampledNOData_[i].resize(chi_param_.size());
      for (int j=0; j<sampledNOData_[i].size();j++)
      {
         sampledNOData_[i][j].resize(Zeta_param_.size());
         for (int k=0; k<sampledNOData_[i][j].size();k++)
         {
            sampledNOData_[i][j][k].resize(Z_param_.size());
         }
      }
   }
   
   //add by jiang .定义混合分数方差、混合分数及进程变量维度的表的存储list
   sampledPVMajorData_.resize(tableNames_.size());      //sampledPVData的第一个维度为cantera表的维度（变量名的数量即温度，化学物质名，混合分数）。

   for (int i=0; i<sampledPVMajorData_.size(); i++)
   {
      sampledPVMajorData_[i].resize(Zeta_param_.size());     //变量的第二个维度定义为混合分数的标准差
      for (int j=0; j<sampledPVMajorData_[i].size();j++)
      {
         sampledPVMajorData_[i][j].resize(Z_param_.size());    //变量的第三个维度定义为混合分数
         for (int k=0; k<sampledPVMajorData_[i][j].size();k++)
         {
            sampledPVMajorData_[i][j][k].resize(PV_Major_param_.size());     //变量的第四个维度定义为进程变量
         }
      }
   }                   //建立三个索引，第一维混合分数的标准差，第二维混合分数，第三维进程变量
	
	//add by wzy .定义混合分数方差、混合分数及进程变量NO维度的表的存储list
   sampledPVNOData_.resize(tableNames_.size());      //sampledPVData的第一个维度为cantera表的维度（变量名的数量即温度，化学物质名，混合分数）。

   for (int i=0; i<sampledPVNOData_.size(); i++)
   {
      sampledPVNOData_[i].resize(Zeta_param_.size());     //变量的第二个维度定义为混合分数的标准差
      for (int j=0; j<sampledPVNOData_[i].size();j++)
      {
         sampledPVNOData_[i][j].resize(Z_param_.size());    //变量的第三个维度定义为混合分数
         for (int k=0; k<sampledPVNOData_[i][j].size();k++)
         {
            sampledPVNOData_[i][j][k].resize(PV_NO_param_.size());     //变量的第四个维度定义为进程变量
         }
      }
   }                   //建立三个索引，第一维混合分数的标准差，第二维混合分数，第三维进程变量
	
   // find the correct FileName to read
   for (int numChi=0; numChi<chi_param_.size();numChi++)
   {
      std::ostringstream chiName;
      chiName << chi_param_[numChi];

      fileName currentFile = canteraFileName + "_" + chiName.str() + ".csv";
      Info << "Reading: " << currentFile << endl;

      // Read tables
      read(currentFile);

      // flip the table if necessary
      List<scalar> hY_(singleData_[tableNames_["Z"]].size(), 0.0);

      if (singleData_[tableNames_["Z"]][singleData_[tableNames_["Z"]].size()-1] <= 0.1)
      {
         for (int i=0; i<singleData_.size(); i++)
         {
            hY_ = 0;
            for (int j=0;j<singleData_[tableNames_["Z"]].size();j++)
            {
               hY_[j] = singleData_[i][singleData_[tableNames_["Z"]].size() - 1 - j];
            }
            singleData_[i] = hY_;
         }
      }

      // Special treatment for boundaries
      singleData_[tableNames_["Z"]][0] = 0;
      singleData_[tableNames_["Z"]][singleData_[tableNames_["Z"]].size()-1] = 1;

      // Calculate hs
      calculateEnthalpy();

      // Calculate Z
      calculateZ();

      for (int numZeta=0; numZeta<Zeta_param_.size();numZeta++)
      {
         // Write Integrated Data in sampledData
         integratedData_ = singleData_;
         //Info << "output singleData" << endl;
         betaPDFIntegration(numChi, Zeta_param_[numZeta]);
         //Info << "betaPDFIntegration" << endl;
         interpolateData();
         //Info << "interpolateData" << endl;
         for (int k=0; k<integratedData_.size();k++)
         {
            sampledNOData_[k][numChi][numZeta] = integratedData_[k];
            sampledMajorData_[k][numChi][numZeta] = integratedData_[k]; //added by wzy把进行pdf积分处理后的标量存到4维表里，用于后面的换维过程
         }
         Info << "output integratedData" << endl;
      }
   }
   Info << "output integratedData  end and success" << endl;
   interpolatePVMajorData();   //调用换维函数，将标量耗散率换成进程变量。
   Info << "interpolatePVMajorData success" << endl;
   interpolatePVNOData();   //调用换维函数，将标量耗散率换成进程变量。

  //Info << "output whole Data" << endl;
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::canteraReader::~canteraReader()
{}

// * * * * * * * * * * * * * * * * MemberFunctions  * * * * * * * * * * * * * * * //

Foam::hashedWordList Foam::canteraReader::getNames()
{
	return tableNames_;
}


void	Foam::canteraReader::write(const int& i,Foam::IOdictionary& dictionary, Foam::OFstream& output)
{
	word dictionaryName=tableNames_[i]+"_table";
	List<List<List<scalar> > >lists=sampledPVMajorData_[i];       //add by jiang，将FPV表分成各个独立的表。
	dictionary.set(dictionaryName,lists);
	dictionary.writeHeader(output);
	output<<dictionaryName<<lists<<";";
}

void	Foam::canteraReader::forewrite(const int& i,Foam::IOdictionary& dictionary, Foam::OFstream& output)
{
	word dictionaryName=tableNames_[i]+"_foretable";
	List<List<List<scalar> > >lists=sampledMajorData_[i];       //add by jiang，将FPV表分成各个独立的表。
	dictionary.set(dictionaryName,lists);
	dictionary.writeHeader(output);
	output<<dictionaryName<<lists<<";";
}

void	Foam::canteraReader::writeNO(const int& i,Foam::IOdictionary& dictionary, Foam::OFstream& output)
{
	word dictionaryName=tableNames_[i]+"_table";
	List<List<List<scalar> > >lists=sampledPVNOData_[i];       //added by wzy，将NO表分成各个独立的表。
	dictionary.set(dictionaryName,lists);
	dictionary.writeHeader(output);
	output<<dictionaryName<<lists<<";";
}
void	Foam::canteraReader::forewriteNO(const int& i,Foam::IOdictionary& dictionary, Foam::OFstream& output)
{
	word dictionaryName=tableNames_[i]+"_foretable";
	List<List<List<scalar> > >lists=sampledNOData_[i];       //add by jiang，将FPV表分成各个独立的表。
	dictionary.set(dictionaryName,lists);
	dictionary.writeHeader(output);
	output<<dictionaryName<<lists<<";";
}

