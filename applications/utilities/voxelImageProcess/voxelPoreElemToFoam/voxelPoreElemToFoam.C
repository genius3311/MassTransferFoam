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
//!   map images, used to initialize Openfoam fields


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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
	
    argList::validArgs.append("voxelPoreElemHeader");
    //argList::validArgs.append("alphaHeader");
    //argList::validOptions.insert("invertAlpha","");
    //argList::validOptions.insert("nGrowAlpha", "label");

	
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"

  
    word voxelPoreElemHeader(args.argRead<word>(1));
    Info<<endl<<"voxelPoreElemHeader:    "<<voxelPoreElemHeader<<endl;

	//voxelImage vximage(voxelPoreElemHeader); 
	voxelImageT<signed int> vximage(voxelPoreElemHeader);

	int3 n=vximage.size3();
	dbl3 xmin=vximage.X0(); //xmin*=1e-6;
	dbl3 dx=vximage.dx();   //dx*=1e-6;

	if (!n[0]) Info<<"\nError: no image read\n"<<endl;
	//if (args.optionFound("invertAlpha"))	vximage.threshold101(1,255);

	runTime.setTime(timeDirs[timeDirs.size()-1], 0);


	const fvBoundaryMesh& boundary = mesh.boundary();
	
	dimensionedScalar PEN ("PEN", dimless, 0.0);
	volScalarField poreElem
	(
	    IOobject
	    (
		"poreElem",
		runTime.timeName(),
		mesh,
		//IOobject::MUST_READ,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh,
	    PEN
	);	
		

	const vectorField & C = mesh.C().internalField();
	//double sumAlpha=0., sumWalpha=1e-64;
	if (vximage.nx())
	{
	  forAll(C,c)
	  {
		int i=(C[c][0]-xmin[0])/dx[0]*0.999999999999;
		int j=(C[c][1]-xmin[1])/dx[1]*0.999999999999;
		int k=(C[c][2]-xmin[2])/dx[2]*0.999999999999;
		poreElem[c]=vximage(i,j,k);
	  }
	}

	forAll(boundary, bi)
		poreElem.boundaryFieldRef()[bi]==poreElem.boundaryField()[bi].patchInternalField();

    //poreElem=fvc::average(linearInterpolate(poreElem));
    //forAll(boundary, bi)	poreElem.boundaryFieldRef()[bi]==poreElem.boundaryField()[bi].patchInternalField();
    
    //while (min(poreElem).value()<2.0)
    /*{
	    for (label itr=0;itr<1;++itr)
	    {
		poreElem=fvc::average(linearInterpolate(poreElem));
		forAll(boundary, bi)	poreElem.boundaryFieldRef()[bi]==poreElem.boundaryField()[bi].patchInternalField();				
	    }
	    poreElem.correctBoundaryConditions();
	    // 对内部场进行四舍五入
	    forAll(poreElem.internalField(), celli) 
	    {
		poreElem[celli] = std::round(poreElem[celli]);
	    }
	    // 更新边界条件
	    forAll(poreElem.boundaryField(), patchi) 
	    {
		fvPatchField<scalar>& pf = poreElem.boundaryFieldRef()[patchi];
		forAll(pf, facei) {
		    pf[facei] = std::round(pf[facei]);
		}
		pf == pf.patchInternalField(); // 保持边界一致性
	    }
    }*/
    //Info<< xmin(poreElem) << "end" << endl;
    
    volScalarField unassigned
	(
	    IOobject
	    (
		"unassigned",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("0", dimless, 0.0)
	);

	// 标记未赋值单元
	forAll(poreElem, c)
	{
	    if (poreElem[c] < SMALL)  // SMALL=1e-15，避免浮点误差
	    {
		unassigned[c] = 1.0; // 标记为需要填充
	    }
	}
	// 获取网格的cell-cell邻接表
	const labelListList& cellCells = mesh.cellCells();

	// 填充未赋值单元
	forAll(unassigned, c)
	{
	    if (unassigned[c] > 0.5)
	    {
		const labelList& neighbours = cellCells[c];
		// 遍历所有相邻单元
		forAll(neighbours, nbrIdx)
		{
		    const label nbrCell = neighbours[nbrIdx];
		    if (poreElem[nbrCell] > SMALL)
		    {
		        poreElem[c] = poreElem[nbrCell];
		    }
		}
	    }
	}

    poreElem.correctBoundaryConditions();
    
    OFstream POF(poreElem.time().timeName()+"/poreElem");
    poreElem.writeHeader(POF);
    poreElem.writeData(POF);
    
    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
