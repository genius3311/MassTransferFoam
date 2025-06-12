/*-------------------------------------------------------------------------*\
 Copyright (C) 2010-2020  Ali Qaseminejad Raeini 

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/

//! Description:
//!   This app converts an unstructured openfoam mesh into a 3D image.


#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"

#include "OFstream.H"

#include "voxelImage.h"
 
using namespace Foam;



int main(int argc, char *argv[])
{
    // 添加命令行选项
    argList::validOptions.insert("alpha", "fieldName");
    argList::validOptions.insert("time", "timeValue | first | last");
    
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"


    IOdictionary meshingDict
    (
        IOobject
        (
            "meshingDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word headerName(meshingDict.lookupOrDefault("headerName",word("")));
    word outputFormat(meshingDict.lookupOrDefault("outputFormat",word("binary")));
    
    // 在读取 alpha 场之前，确保用户了解当前选择的时间步
    Info << "Reading data from time: " << runTime.timeName() << nl;
    if (args.optionFound("time"))
    {
        Info << "Time was explicitly set by -time option" << nl;
    }
    else
    {
        Info << "No -time option provided, using latest time step by default" << nl;
    }
    // 处理 alpha 场选项
    word alphaName = "alpha1";  // 默认场名称
    if (args.optionFound("alpha"))
    {
        args.optionLookup("alpha")() >> alphaName;
        Info << "Using alpha field: " << alphaName << nl;
    }
    else
    {
        Info << "Alpha field name not specified, using default: alpha1" << nl;
    }

    // 检查 alpha 场是否存在
    bool useAlpha = false;
    //volScalarField alphaField;
    tmp<volScalarField> talphaField;
    
    IOobject alphaHeader
    (
        alphaName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    
    if (alphaHeader.typeHeaderOk<volScalarField>(true))
    {
        //alphaField = volScalarField(alphaHeader, mesh);
        talphaField = tmp<volScalarField>(new volScalarField(alphaHeader, mesh));
        useAlpha = true;
        Info << "Alpha field found: " << alphaName << ". Will be used to distinguish phases." << nl;
    }
    else
    {
        Info << "Warning: Alpha field not found: " << alphaName << ". Will generate single-phase pore space." << nl;
    }

    int3 nnn;
    dbl3 xmin, dx;
    voxelImage rock(headerName);
    surfaceVectorField Cf =  mesh.Cf();
    surfaceVectorField Sf =  mesh.Sf();

    if (rock.nx()==0)
    {
		if(headerName.size()) Info<<"\n\nWarning: ignoring header file from system/meshingDict:"<<headerName<<endl<<endl;
		headerName="";
		Info <<"computing mesh extents"<<endl;
		dx[0]=std::pow(gAverage(mesh.V()),1./3.); dx[1]=dx[0]; dx[2]=dx[1];
		rock.dxCh()=dx;

		boundBox box(mesh.points());
		rock.X0Ch()[0]=box.min()[0];  rock.X0Ch()[1]=box.min()[1]; rock.X0Ch()[2]=box.min()[2];
		//rock.X0Ch()[1]-=dx[1]; rock.X0Ch()[2]-=dx[2];
		nnn[0]=((box.max()[0]-box.min()[0]) + 0.6*dx[0])/dx[0];
		nnn[1]=((box.max()[1]-box.min()[1]) + 0.6*dx[1])/dx[1];
		nnn[2]=((box.max()[2]-box.min()[2]) + 0.6*dx[2])/dx[2];

		Info<<"X0: "<<rock.X0()[0]<<" "<<rock.X0()[1]<<" "<<rock.X0()[2]<<endl;
		Info<<"dx:   "<<rock.dx()[0]<<"  "<<rock.dx()[1]<<"  "<<rock.dx()[2]<<endl;
		Info<<"nnn:    "<<nnn[0]<<" "<<nnn[1]<<" "<<nnn[2]<<endl;


		const dbl3 X0 = rock.X0();
		const dbl3 dx = rock.dx();
		label patchI = mesh.boundaryMesh().findPatchID("Left");
		if (patchI >= 0)
		/*{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[0];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[0]>.99 && mag(avgPos-X0[0]) < 2.*dx[0])
				rock.X0Ch()[0] = avgPos;
		}*/
		{
                       double avgPos = gAverage(Cf.boundaryField()[patchI])[0];
                       vector avgN = gSum(Cf.boundaryField()[patchI]); 
                       scalar magN = mag(avgN);
                       if (magN > SMALL) // 添加安全检查，避免除以零
                       {
                             avgN /= magN;
      			      if (avgN[0]>.99 && mag(avgPos-X0[0]) < 2.*dx[0])
          		      rock.X0Ch()[0] = avgPos;
 			}
 			else
    			{
    			    Info << "Warning: Zero magnitude vector for patch " << mesh.boundaryMesh()[patchI].name() << endl;
   			 }
		}

		patchI = mesh.boundaryMesh().findPatchID("Right");
		if (patchI >= 0)
		/*{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[0];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[0]>.99 && mag(avgPos-X0[0]-nnn[0]*dx[0]) < 2.*dx[0])
				rock.dxCh()[0] = (avgPos-X0[0])/nnn[0];
		}*/
		{
                       double avgPos = gAverage(Cf.boundaryField()[patchI])[0];
                       vector avgN = gSum(Cf.boundaryField()[patchI]); 
                       scalar magN = mag(avgN);
                       if (magN > SMALL) // 添加安全检查，避免除以零
                       {
                             avgN /= magN;
      			      if (avgN[0]>.99 && mag(avgPos-X0[0]-nnn[0]*dx[0]) < 2.*dx[0])
          		      rock.dxCh()[0] = (avgPos-X0[0])/nnn[0];
 			}
 			else
    			{
    			    Info << "Warning: Zero magnitude vector for patch " << mesh.boundaryMesh()[patchI].name() << endl;
   			 }
		}

		patchI = mesh.boundaryMesh().findPatchID("Bottom");
		if (patchI >= 1)
		/*{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[1];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[1]>.99 && mag(avgPos-X0[1]) < 2.*dx[1])
				rock.X0Ch()[1] = avgPos;
		}*/
		{
                       double avgPos = gAverage(Cf.boundaryField()[patchI])[1];
                       vector avgN = gSum(Cf.boundaryField()[patchI]); 
                       scalar magN = mag(avgN);
                       if (magN > SMALL) // 添加安全检查，避免除以零
                       {
                             avgN /= magN;
      			      if (avgN[1]>.99 && mag(avgPos-X0[1]) < 2.*dx[1])
          		      rock.X0Ch()[1] = avgPos;
 			}
 			else
    			{
    			    Info << "Warning: Zero magnitude vector for patch " << mesh.boundaryMesh()[patchI].name() << endl;
   			 }
		}

		patchI = mesh.boundaryMesh().findPatchID("Top");
		if (patchI >= 1)
		/*{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[1];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[1]>.99 && mag(avgPos-X0[1]-nnn[1]*dx[1]) < 2.*dx[1])
				rock.dxCh()[1] = (avgPos-X0[1])/nnn[1];
		}*/
		{
                       double avgPos = gAverage(Cf.boundaryField()[patchI])[1];
                       vector avgN = gSum(Cf.boundaryField()[patchI]); 
                       scalar magN = mag(avgN);
                       if (magN > SMALL) // 添加安全检查，避免除以零
                       {
                             avgN /= magN;
      			      if (avgN[1]>.99 && mag(avgPos-X0[1]) < 2.*dx[1])
          		      rock.dxCh()[1] = (avgPos-X0[1])/nnn[1];
 			}
 			else
    			{
    			    Info << "Warning: Zero magnitude vector for patch " << mesh.boundaryMesh()[patchI].name() << endl;
   			 }
		}

		patchI = mesh.boundaryMesh().findPatchID("Back");
		if (patchI >= 2)
		/*{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[2];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[2]>.99 && mag(avgPos-X0[2]) < 2.*dx[2])
				rock.X0Ch()[2] = avgPos;
		}*/
		{
                       double avgPos = gAverage(Cf.boundaryField()[patchI])[2];
                       vector avgN = gSum(Cf.boundaryField()[patchI]); 
                       scalar magN = mag(avgN);
                       if (magN > SMALL) // 添加安全检查，避免除以零
                       {
                             avgN /= magN;
      			      if (avgN[2]>.99 && mag(avgPos-X0[2]) < 2.*dx[2])
          		      rock.X0Ch()[2] = avgPos;
 			}
 			else
    			{
    			    Info << "Warning: Zero magnitude vector for patch " << mesh.boundaryMesh()[patchI].name() << endl;
   			 }
		}
		
		patchI = mesh.boundaryMesh().findPatchID("Front");
		if (patchI >= 0)
		/*{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[2];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[2]>.99 && mag(avgPos-X0[2]-nnn[2]*dx[2]) < 2.*dx[2])
				rock.dxCh()[2] = (avgPos-X0[2])/nnn[2];
		}*/
		{
                       double avgPos = gAverage(Cf.boundaryField()[patchI])[2];
                       vector avgN = gSum(Cf.boundaryField()[patchI]); 
                       scalar magN = mag(avgN);
                       if (magN > SMALL) // 添加安全检查，避免除以零
                       {
                             avgN /= magN;
      			      if (avgN[2]>.99 && mag(avgPos-X0[2]-nnn[2]*dx[2]) < 2.*dx[2])
          		      rock.dxCh()[2] = (avgPos-X0[2])/nnn[2];
 			}
 			else
    			{
    			    Info << "Warning: Zero magnitude vector for patch " << mesh.boundaryMesh()[patchI].name() << endl;
   			}
		}


		Info<<"->X0:   "<<rock.X0()[0]<<"  "<<rock.X0()[1]<<"  "<<rock.X0()[2]<<endl;
		Info<<"->dx:   "<<rock.dx()[0]<<"  "<<rock.dx()[1]<<"  "<<rock.dx()[2]<<endl;
		double dxAvg = (rock.dxCh()[0]+rock.dxCh()[1]+rock.dxCh()[2])/3.;

		rock.dxCh()[0]=dxAvg;  rock.dxCh()[1]=dxAvg;  rock.dxCh()[2]=dxAvg;
		Info<<"->dx:   "<<rock.dx()[0]<<"  "<<rock.dx()[1]<<"  "<<rock.dx()[2]<<endl;

		//rock.X0Ch()*=1e6; 
		//rock.dxCh()*=1e6; 
	}   
	else
       {
		Info <<"mesh size read from header"<<endl;
		nnn=rock.size3();
		rock.reset(0,0,0,0);

	}


	xmin=rock.X0(); dx=rock.dx();
	//dx*=1e-6; 
	//xmin*=1e-6; 

	Info<<"xmin: "<<xmin[0]<<" "<<xmin[1]<<" "<<xmin[2]<<endl;
	Info<<"dx:   "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;
	Info<<"nnn:    "<<nnn[0]<<" "<<nnn[1]<<" "<<nnn[2]<<endl;


	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);


	//voxelField<double> pVoxel(nnn[0],nnn[1],nnn[2],0.);
	//rock.reset(nnn[0],nnn[1],nnn[2],1);
	// 初始化整个体素场为骨架值 (0)
        rock.reset(nnn[0], nnn[1], nnn[2], 0);
        voxelImage rockAlpha = rock;

	const vectorField & C =	mesh.C();
	const scalarField & V =	mesh.V();

	Info<<"mesh size: "<<C.size()<<endl;

       /*forAll(C,c)
	{
		int i=(C[c][0]-xmin[0])/dx[0];
		int j=(C[c][1]-xmin[1])/dx[1];
		int k=(C[c][2]-xmin[2])/dx[2];

		if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
		else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
	}*/
       forAll(C, c)
       {
           int i = (C[c][0]-xmin[0])/dx[0];
           int j = (C[c][1]-xmin[1])/dx[1];
           int k = (C[c][2]-xmin[2])/dx[2];

           if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])
           {
               rock(i,j,k)=1;
               // 使用 alpha 值来区分两种孔隙相
               const volScalarField& alphaField = talphaField();
               //rockAlpha(i,j,k) = alphaField[c];
               if (useAlpha)
               {                  
                   if (alphaField[c] > 0.5)          
                       rockAlpha(i,j,k) = 1;  // alpha ≈ 1 的区域，设为 1
                   else
                       rockAlpha(i,j,k) = 2;  // alpha ≈ 0 的区域，设为 2
               }
               else
                   rockAlpha(i,j,k) = 1;  // 没有 alpha 场时，所有孔隙统一用 1 表示
               
           }
           else
               Info << "Error: ijk: " << i << " " << j << " " << k << endl;
       }
       // 添加额外的填充点以确保所有网格单元都被捕获
       forAll(C,c)
       {
           double dxi = std::pow(V[c],1./3.);
           // 循环遍历8个角点
	    for (int cx = -1; cx <= 1; cx += 2)
	    {
		for (int cy = -1; cy <= 1; cy += 2)
		{
		    for (int cz = -1; cz <= 1; cz += 2)
		    {
		        int i = (C[c][0] - xmin[0] + cx * dxi/4.) / dx[0];
		        int j = (C[c][1] - xmin[1] + cy * dxi/4.) / dx[1];
		        int k = (C[c][2] - xmin[2] + cz * dxi/4.) / dx[2];
		        
		        if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])
		        {
		            rock(i,j,k)=1;
		            const volScalarField& alphaField = talphaField();
		            //rockAlpha(i,j,k) = alphaField[c];
		            if (useAlpha)
			    {
				
				if (alphaField[c] > 0.5)             
				    rockAlpha(i,j,k) = 1;  // alpha ≈ 1 的区域，设为 1
				else
				    rockAlpha(i,j,k) = 2;  // alpha ≈ 0 的区域，设为 2
			    }
			    else
			    {
				rockAlpha(i,j,k) = 1;  // 没有 alpha 场时，所有孔隙统一用 1 表示
			    }
		        }
		    }
		}
	    }
       }
       
       /*forAll(C,c)
       {
           double dxi = std::pow(V[c],1./3.);
           // 循环遍历8个角点
	    for (int cx = -1; cx <= 1; cx += 2)
	    {
		for (int cy = -1; cy <= 1; cy += 2)
		{
		    for (int cz = -1; cz <= 1; cz += 2)
		    {
		        int i = (C[c][0] - xmin[0] + cx * dxi/4.) / dx[0];
		        int j = (C[c][1] - xmin[1] + cy * dxi/4.) / dx[1];
		        int k = (C[c][2] - xmin[2] + cz * dxi/4.) / dx[2];
		        
		        if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])
		            rock(i,j,k)=1;
		    }
		}
	    }
       }*/
	
	
	Info << "Applying modified smoothing to preserve phase interfaces..." << nl;
	rock.FaceMedian06(2,4);
	forAll(Cf,c)
	{
		int i=(Cf[c][0]-xmin[0])/dx[0];
		int j=(Cf[c][1]-xmin[1])/dx[1];
		int k=(Cf[c][2]-xmin[2])/dx[2];

		if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=1;
	}
	rock.FaceMedian06(1,5);
	rock.FaceMedian06(1,5);
	rock.FaceMedian06(1,5);
	
	//voxelImage rockAlpha = rock;
	//rockAlpha = rock;
	/*forAll(C, c)
        {
           int i = (C[c][0]-xmin[0])/dx[0];
           int j = (C[c][1]-xmin[1])/dx[1];
           int k = (C[c][2]-xmin[2])/dx[2];
           
           // 使用 alpha 值来区分两种孔隙相
           const volScalarField& alphaField = talphaField();

           if (rock(i,j,k) > 0)
           //if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])
           {
               //if (rock(i,j,k) > 0)
               if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])
               {                                    
                   if (alphaField[c] > 0.5)          
                       rockAlpha(i,j,k) = 1;  // alpha ≈ 1 的区域，设为 1
                   else
                       rockAlpha(i,j,k) = 2;  // alpha ≈ 0 的区域，设为 2
                   //Info<< "TTTTTTTTTTTTTTTest" << endl;
               }
               else
                   rockAlpha(i,j,k) = 0;
           }
           else
               Info << "Error: ijk: " << i << " " << j << " " << k << endl;
        }
        */
       
       // check the interface again
       for (label i=0; i<nnn[0]; i++)
        {
            for (label j=0; j<nnn[1]; j++)
            {
                for (label k=0; k<nnn[2]; k++)
                {
                    if (rock(i,j,k) < 1 && rockAlpha(i,j,k)>0)
                        rockAlpha(i,j,k) = rock(i,j,k);
                    if (rockAlpha(i,j,k) < 1 && rock(i,j,k)>0)
                    {
                        double sum = 0.0;
                        int count = 0;
                        for (int ni = -1; ni <= 1; ni++) {
                            for (int nj = -1; nj <= 1; nj++) {
                                for (int nk = -1; nk <= 1; nk++) {
                                     if (ni == 0 && nj == 0 && nk == 0)     continue;
                                     if (rock(i+ni, j+nj, k+nk) > 0) {
                                         sum += rock(i+ni, j+nj, k+nk);
                                         count++;
                                     }
                                }
                            }
                        }
			 if (count > 0)   rockAlpha(i,j,k) = sum / count;
                    }
                }
            }
        }


	rock.replacexLayer(0,1);
	rock.replacexLayer(nnn[0]-1,nnn[0]-2);

	rock.replaceyLayer(0,1);
	rock.replaceyLayer(nnn[1]-1,nnn[1]-2);

	rock.replacezLayer(0,1);
	rock.replacezLayer(nnn[2]-1,nnn[2]-2);
	
	
	rockAlpha.replacexLayer(0,1);
	rockAlpha.replacexLayer(nnn[0]-1,nnn[0]-2);

	rockAlpha.replaceyLayer(0,1);
	rockAlpha.replaceyLayer(nnn[1]-1,nnn[1]-2);

	rockAlpha.replacezLayer(0,1);
	rockAlpha.replacezLayer(nnn[2]-1,nnn[2]-2);

	if (outputFormat=="binary")
	{
		 rock.write("vxlImage"+imgExt());
	}
	else
	{
		 rock.write("vxlImage.dat");
        }
        // 输出体素数据
        if (outputFormat=="binary")
        {
            //rockAlpha.write(runTime.path()+"/vxlImageAlpha"+imgExt());
            rockAlpha.write(runTime.timeName()+"-vxlImageAlpha"+imgExt());         
        }
        else
        {
            rockAlpha.write(runTime.timeName()+"-vxlImageAlpha.dat");
        }
     


    Info<< "\nend" << endl;

    return 0;
}


