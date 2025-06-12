/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvmLaplacian.H"
#include "fvCFD.H"
#include "unitConversion.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf,
    volVectorField::Boundary& gradAlphab
    //volScalarField::Boundary& alphaSb
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                //boundary[patchi].nf()
                nw_.boundaryField()[patchi]
                
            );
            
            //method 1
            /*vectorField gAlphaS = gradAlphab[patchi];
            primitivePatchInterpolation pinterpolator(mesh.boundaryMesh()[patchi]);
            gAlphaS = 0.1*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(gAlphaS) ) + 0.9*gradAlphab[patchi];
            vectorField nsp = gAlphaS/(mag(gAlphaS) + deltaN_.value());
            
            vectorField ss(nsp - (nsp&nf)*nf);
            ss/=mag(ss)+1e-37;
            nsp=sin(theta)*ss+cos(theta)*nf;
            nHatp = pinterpolator.faceToPointInterpolate(nsp);
            //nHatp /= (mag(nHatp) + deltaN_.value());
            //nHatp.setOriented();
            //nHatp[patchi] = nsp/(mag(nsp) + deltaN_.value());
            //acap.gradient() == 0.5*(1.-acap*acap)*((boundary[patchi].nf() & nsp)*mag(gradAlphab[patchi]));
            acap.gradient() = (boundary[patchi].nf() & nHatp)*mag(gradAlphab[patchi]);
            acap.evaluate();*/
            /*alphaContactAngleTwoPhaseFvPatchScalarField& alphaSbCap =
		const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
		( 
		    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
		    (
		        alphaSb[patchi]
		    ) 
		);            
	    alphaSbCap.gradient() = (boundary[patchi].nf() & nsp)*mag(gAlphaS);
	    alphaSbCap.gradient() == (boundary[patchi].nf() & nsp)*mag(gAlphaS);
	    gradAlphab[patchi] = nsp*mag(gradAlphab[patchi]);
	    alphaSbCap.evaluate();*/
	    
	    //method 2
	    /*vectorField ss(nHatp - (nHatp&nf)*nf);
	    ss/=mag(ss)+1e-37;
	    nHatp=sin(theta)*ss+cos(theta)*nf;
	    nHatp /= (mag(nHatp) + deltaN_.value());*/

            // the original calculation method
            // Reset nHatp to correspond to the contact angle            
            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    //init alpha smoothed
    //volScalarField alpha1s = alpha1_;
    volScalarField alpha1s = ( min( max(alpha1_*1.02-0.01,1e-3), (1.-(1e-3)) ) );

    // smooth alpha1 by successive interpolation to face and cell
    /*volScalarField alpha1Stmp = alpha1s;
    surfaceScalarField alpha1f_ = linearInterpolate(min( max(alpha1_,0.), 1. ));
    volScalarField a1a2 = 0.95+0.1*sqrt(alpha1Stmp*(1.-alpha1Stmp));
    alpha1Stmp = a1a2*alpha1Stmp+(1.-a1a2)*fvc::average(alpha1f_);
    alpha1Stmp.correctBoundaryConditions();*/
    /*a1a2 = 2.*sqrt(alpha1Stmp*(1.-alpha1Stmp));
    alpha1Stmp = a1a2*alpha1Stmp + (1.-a1a2)*fvc::average(fvc::interpolate(alpha1Stmp)); 
    alpha1Stmp.correctBoundaryConditions();
    alpha1Stmp = 0.99*alpha1Stmp + 0.01*fvc::average(fvc::interpolate(alpha1Stmp));
    //alpha1s = atanh(1.8*alpha1Stmp-0.9); //- SYN123433
    */
    //alpha1s = alpha1Stmp;
    // Raeini's thesis, equation (2.14)
    for (int i=0;i<nSK_;i++)
    {
        alpha1s = cSK_ * fvc::average(fvc::interpolate(alpha1s)) + (1.0 - cSK_) * alpha1s;
    }
    alpha1s.correctBoundaryConditions();
    volScalarField a1a2 = 2.*sqrt(alpha1s*(1.-alpha1s));
    
    // Cell gradient of alpha
    //const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));  
    //const volVectorField gradAlpha(fvc::grad(alpha1s, "nHat"));
    //volVectorField gradAlpha(fvc::grad(alpha1s, "nHat"));
    volVectorField gradAlpha = fvc::reconstruct(fvc::snGrad(alpha1s)*mesh.magSf());   // used in porefoam
    //gradAlpha.correctBoundaryConditions();
    
    //update interface vector at cell center, used for filtering
    volScalarField magGradAlpha = mag(gradAlpha) + deltaN_;
    //nI_ = gradAlpha/(Foam::mag(gradAlpha) + deltaN_);
    nI_ = a1a2*gradAlpha/magGradAlpha;
    nI_ = nI_/(mag(nI_) + 1e-18);
    /*----------------------------------------------*/

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    //surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    surfaceVectorField nHatfv(fvc::interpolate(nI_,"smoothScheme"));
    nHatfv/=mag(nHatfv)+1e-12;
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    //correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());
    //correctContactAngle(nHatfv.boundaryFieldRef(), gradAlpha.boundaryFieldRef());
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField(), gradAlpha.boundaryFieldRef());
    //correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField(), gradAlpha.boundaryFieldRef(), alpha1s.boundaryFieldRef());
    
    volScalarField a1a2Relaxed = 0.95*a1a2;    ///smoothingRelaxFactor_*a1a2;
    //volScalarField 
    a1xa2 = a1a2*(1.-0.000001)+0.000001;
    for (int i=0; i<nSK_;i++)     //(int i=0; i<smoothingKernel_;i++) 
    {	
	nI_=(1.-a1a2Relaxed)*nI_+a1a2Relaxed*fvc::average(fvc::interpolate(a1a2Relaxed*nI_,"smoothScheme"));	
	nI_.correctBoundaryConditions();
	nI_=(1.-a1a2Relaxed)*nI_+a1a2Relaxed*fvc::average(fvc::interpolate(nI_));	
	nI_.correctBoundaryConditions();
	nHatfv = fvc::interpolate(nI_,"smoothScheme");
	nHatfv = nHatfv/(mag(nHatfv) + 1e-12);
	//correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());
	correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField(), gradAlpha.boundaryFieldRef());	
	//correctContactAngle(nHatfv.boundaryFieldRef(),gradAlpha.boundaryFieldRef(),nS_.boundaryFieldRef(),alpha1S_.boundaryFieldRef());
    }
    nI_ = nI_/(mag(nI_) + 1e-12);
    

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;
    
    surfaceScalarField snGradAlpha = fvc::snGrad(alpha1s);
    forAll(nHatf_,i)  if (nHatf_[i]*snGradAlpha[i] < 0) nHatf_[i]*=0.9;
    
    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    //K_ = -fvc::div(nHatf_) + (nI_ & fvc::grad(nHatfv) & nI_);

    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
    
    nHatfv = fvc::interpolate(nI_);  // the scheme is different from above
    nHatf_ = nHatfv & Sf;
    forAll(nHatf_,i)   if (nHatf_[i]*snGradAlpha[i] < 0) nHatf_[i]*=-0.1;
    volScalarField alphash = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);
    surfaceScalarField delS= fvc::snGrad(alphash);
    if (nSK_) ///. smoothing Kc
    {
	    // this correction helps stablizing also improves the accuracy for capillary pressure,
	    // it has some theoretical basis but the coefficients here are chosen emperically
	    volScalarField mgK=mag(K_);
	    K_/=mag( 1. + 0.1*((alpha1s-0.5)+0.1*K_/(magGradAlpha+5.*mgK)) * K_/(magGradAlpha+2.*mgK) )+1e-12;
	    K_.correctBoundaryConditions();
	    surfaceScalarField WKf=0.02+0.08*fvc::interpolate(a1xa2)+(mag(delS)/(mag(delS)+deltaN_));
	    K_ = fvc::average(fvc::interpolate(K_)*WKf )/fvc::average(WKf);
	    K_ = fvc::average(fvc::interpolate(K_)*WKf )/fvc::average(WKf);
	    K_.correctBoundaryConditions();
    }   
}

//Re-calculate capillary flux
void Foam::interfaceProperties::calculatePhic()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceScalarField& magSf = mesh.magSf();
    const fvBoundaryMesh& boundary = mesh.boundary();
 
    volScalarField alpha_pc
    (
        IOobject 
        (
            "alpha_pc",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
     alpha1_,
     pc_.boundaryField().types()
    );

    // Sharpen interface function
    // Raeini's thesis, equation (2.22) and (2.23)
    alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);

    alpha_pc.correctBoundaryConditions();

    surfaceScalarField deltasf = fvc::snGrad(alpha_pc);

    //surface tension force
    //surfaceScalarField stf = fvc::interpolate(sigmaK())*deltasf;    
    surfaceScalarField stf = fvc::interpolate(sigmaK()*a1xa2)/fvc::interpolate(a1xa2)*deltasf;


    //surface tension force flux
    phic_ = stf*magSf;
    
    // prepare for filter, dont know the reason for "Operator - : undefined for oriented and unoriented types"
    sgPc_ = phic_;
    sgPc_.setOriented();
    phic_.setOriented();
    // Filtering capillary forces parallel to interfaces
    surfaceScalarField dInterf=(mag(deltasf)/(mag(deltasf)+1e-4*deltaN_));
    sgPcErr_ *= dInterf;
    sgPcErr_.setOriented();
    phic_ = sgPc_ - sgPcErr_;
    
    pc_ = 0.9*pc_+0.1*fvc::average(fvc::interpolate(pc_));// hoping to avoid solver crash for bad quality meshes

    for(int nonOrth=0; nonOrth<=nNonOrthogonalCorrectors_; nonOrth++)
    {
		//solve for pc
        fvScalarMatrix pcEqn
        (
            fvm::laplacian( pc_) == fvc::div(phic_)
        );

        pcEqn.setReference(pcRefCell_, pcRefValue_);

        solverPerformance residual = pcEqn.solve();

        
	if (nonOrth==0) eqnResidual_ = residual.initialResidual();

		//add flux of pc to capillary flux
        if (nonOrth == nNonOrthogonalCorrectors_)
        { 
            phic_-=pcEqn.flux();
        }
        pc_.correctBoundaryConditions();
    }
    
    // Filtering capillary forces parallel to interfaces
    // update sgPcErr_ with time
    
    surfaceVectorField nHatfv(fvc::interpolate(nI_,"smoothScheme"));
    nHatfv/=mag(nHatfv)+1e-12;
    
    volVectorField gPc_= fvc::reconstruct((phic_));
    gPc_.correctBoundaryConditions();
    volVectorField vgPcllInterface=(gPc_)-((gPc_)&(nI_))*(nI_);
    surfaceVectorField gPcllInterface=fvc::interpolate(vgPcllInterface); 
    gPcllInterface=(gPcllInterface)-((gPcllInterface)&(nHatfv))*(nHatfv);
    gPcllInterface*=(1.-mag(nHatf_)/magSf);
    sgPcErr_=dInterf*(min(2.*(1.4-mag(nHatf_)/magSf), 1.)*(sgPcErr_+fcCorrectTangent_*(gPcllInterface& mesh.Sf())) );  //fcCorrectTangent_
    surfaceScalarField sigmaf_ = fvc::interpolate(sigmaPtr_->sigma());
    surfaceScalarField stfTypicalLow=mag((sigmaf_*0.1)*deltasf*deltasf*magSf);
    forAll(sgPcErr_,i)
    {
        if (phic_[i]*sgPcErr_[i] < -1e-6*stfTypicalLow[i])  sgPcErr_[i]  *= 0.99;
        if( mag(phic_[i])<0.5*mag(sgPcErr_[i]) ) {   sgPcErr_[i] = 0.5*phic_[i]; }
    }
    forAll(boundary, bI)
    {
        if (!boundary[bI].coupled())
        {
            sgPcErr_.boundaryFieldRef()[bI] = 0.;
        }
        else
        {
            Field<scalar> & psgPcErr_ = sgPcErr_.boundaryFieldRef()[bI];
            Field<scalar> & psgPc_ = sgPc_.boundaryFieldRef()[bI];
            forAll(psgPcErr_,i)
            {
                if(psgPc_[i]*psgPcErr_[i] < -1e-6*stfTypicalLow[i])  psgPcErr_[i] *= 0.99;
                if( mag(psgPc_[i])<0.5*mag(psgPcErr_[i]) ) {  psgPcErr_[i] = 0.5*psgPc_[i];}
            }
        }
    }
    
    
    //Filtering capillary fluxes 
    //Raeini's thesis, equation (2.26)
    //surfaceScalarField sigmaf_ = fvc::interpolate(sigmaPtr_->sigma());
    surfaceScalarField stfThreshold = ((1.*fcdFilter_)*mag(sgPc_)+ mag((sigmaf_*0.1)*deltasf*deltasf*magSf) );
    stfThreshold.setOriented();
    stfThreshold = min(max(phic_, -stfThreshold),stfThreshold);
    phic_-=stfThreshold; 
    
    gPc_= fvc::reconstruct((phic_));
    gPc_.correctBoundaryConditions();
    volVectorField gPcS=fvc::average(fvc::interpolate(gPc_));
    stfThreshold = 2.*mag((fvc::interpolate(gPcS)) & mesh.Sf());
    stfThreshold.setOriented();
    phic_ = min(max(phic_, -stfThreshold),stfThreshold);
    
    gPc_= fvc::reconstruct((phic_));
    gPc_.correctBoundaryConditions();
    phic_ = linearInterpolate(gPc_) & alpha1_.mesh().Sf();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),
    //- smoothing coefficient for alphaS (generally 0.5)
    cSK_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cSK")
        )
    ),
    //-number of smoothing cycle for alphaS
    nSK_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("nSK")
        )
    ),
    //- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
    cPc_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cPc")
        )
    ),
    //-number of smoothing cycle for nw
    wallSmoothingKernel_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("wallSmoothingKernel")
        )
    ),
    fcdFilter_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("fcdFilter") 
        )
    ),
    fcCorrectTangent_
    ( 
        readScalar
        ( 
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("fcCorrectTangent") 
        ) 
    ),
    /*fcCorrectTangentRelax_
    ( 
        readScalar
        ( 
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("fcCorrectTangentRelax") 
        ) 
    ),*/
    
    //-number of non-orthogonal corrector loop
    nNonOrthogonalCorrectors_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("nNonOrthogonalCorrectors")
        )
    ),
    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(alpha1.mesh().V()))
    ),

    alpha1_(alpha1),
    U_(U),
    
    //unit normal vector to the solid wall, smoothed using a Gaussian filter
    nw_
    (
        IOobject
        ( 
            "nw",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),     
        dimensionedVector("nw", dimless, vector(0.,0.,0.))
    ),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    //intreface normal vector at cell center
    nI_
    (
        IOobject
        (
            "nI",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, Zero)
    ),
    //capillary pressure
    pc_
    (
        IOobject
        (
            "pc",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh()
    ),
    gPc_
    (
        IOobject
        ( 
            "gPc",  
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::READ_IF_PRESENT   
        ),
        alpha1_.mesh(),
        dimensionedVector("gPc", dimPressure/dimLength, vector(0.,0.,0.)),
        pc_.boundaryField().types()
    ),
    sgPc_
    (
        IOobject
        ( 
            "sgPc",  
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::READ_IF_PRESENT   
        ), 
        //linearInterpolate(gPc_) & alpha1_.mesh().Sf()
        alpha1_.mesh(),
        dimensionedScalar("sgPc", dimPressure / dimLength*dimArea, 0.0)
    ),
    a1xa2
    (
        IOobject
        (
            "a1xa2",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_
    ),
    
    
    //Specific cell where the capillary pressure value is used as refference
    pcRefCell_
    (
        readLabel
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefCell")
        )
    ),
    //value of the capillary pressure in the specific cell
    pcRefValue_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefValue")
        )
    ),
    
    sgPcErr_
    (  
        IOobject
        ( 
            "sgPce",
            alpha1_.time().timeName(), 
            alpha1_.mesh(), 
            IOobject::READ_IF_PRESENT 
        ), //! NO_WRITE means restarting releases the filters
        alpha1_.mesh(),   
        dimensionedScalar("sgPce", dimPressure/dimLength*dimArea, 0.)
    ),
    
    //capillary flux
    phic_
    (
        IOobject
        (
            "phic",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("phic", dimPressure / dimLength*dimArea, 0.0)
        //linearInterpolate(gPc_) & alpha1_.mesh().Sf()
    )
{
    setRefCell
    (
        pc_,
        alpha1_.mesh().solutionDict().subDict("PIMPLE"),
        pcRefCell_,
        pcRefValue_
    );
    
    const fvMesh& mesh = alpha1_.mesh();
    { // smooth nw
        const fvBoundaryMesh& boundary = mesh.boundary();
        forAll(boundary, bi)
        {
	    if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(alpha1_.boundaryField()[bi]))
	    {
		nw_.boundaryFieldRef()[bi]=boundary[bi].nf();
		nw_.boundaryFieldRef()[bi]==boundary[bi].nf();  //tomakesure
	    }
	}
	forAll(boundary, bi)
	{
	    if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(alpha1_.boundaryField()[bi]))
	    {
		primitivePatchInterpolation pinterpolator(mesh.boundaryMesh()[bi]);
		for (int i=0;i<wallSmoothingKernel_;i++)
		{
			nw_.boundaryFieldRef()[bi]==
			pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(nw_.boundaryField()[bi]));
		}
	    }
	}
	nw_.ref().field() = fvc::interpolate(fvc::average(nw_))->internalField();		
	nw_ /=  mag(nw_) + 1e-15;
    }
    
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::sigmaK() const
{
    return sigmaPtr_->sigma()*K_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::interfaceProperties::correct()
{
    calculateK();
    calculatePhic();
    //#include "separatePc_filter.H"
}


bool Foam::interfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //
