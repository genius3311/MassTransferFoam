/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvPatchFields.H"
#include "heterogeneousAlphaContactAngleFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heterogeneousAlphaContactAngleFvPatchScalarField::
heterogeneousAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    //alphaContactAngleTwoPhaseFvPatchScalarField(p, iF)
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    HCAName_("HCA"),
    theta0_(0.0)
{}

Foam::heterogeneousAlphaContactAngleFvPatchScalarField::
heterogeneousAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    //alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict)
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    HCAName_(dict.getOrDefault<word>("HCA", "HCA")),
    theta0_(dict.get<scalar>("theta0"))
{
    
    evaluate();
}

Foam::heterogeneousAlphaContactAngleFvPatchScalarField::
heterogeneousAlphaContactAngleFvPatchScalarField
(
    const heterogeneousAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    //alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper)
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    HCAName_(gcpsf.HCAName_),
    theta0_(gcpsf.theta0_)
{}



Foam::heterogeneousAlphaContactAngleFvPatchScalarField::
heterogeneousAlphaContactAngleFvPatchScalarField
(
    const heterogeneousAlphaContactAngleFvPatchScalarField& gcpsf
)
:
    //alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf)
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    HCAName_(gcpsf.HCAName_),
    theta0_(gcpsf.theta0_)
{}


Foam::heterogeneousAlphaContactAngleFvPatchScalarField::
heterogeneousAlphaContactAngleFvPatchScalarField
(
    const heterogeneousAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    //alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF)
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    HCAName_(gcpsf.HCAName_),
    theta0_(gcpsf.theta0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::heterogeneousAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    //const fvPatchField<scalar>& rhop = patch().lookupPatchField<const::Foam::volScalarFie ld, scalar>("rho");
    //scalarField phip =  this->patch().lookupPatchField(const::Foam::volScalarFie"rho");
    //scalarField test =pd.boundaryMesh()[patch];

    //const fvPatchField<scalar>& HCA = patch().lookupPatchField<volScalarField, scalar>("HCA");
    //const volScalarField& HCA = this->db().objectRegistry::lookupObject<volScalarField>("HCA");
    //const auto& HCA = patch().lookupPatchField<volScalarField>("HCA");
    /*const volScalarField& HCA = 
    (
      this->db().objectRegistry::template
      lookupObject<volScalarField>("HCA")
     );*/
     
     const volScalarField& alpha1 = this->db().objectRegistry::lookupObject<volScalarField>("alpha.phase1");
     const fvMesh& mesh = alpha1.mesh();
     volScalarField HCA_
     (
        IOobject
        (
            //"HCA",
            HCAName_,
            //mesh.time().timeName(),
            "./0",
            mesh,
            IOobject::MUST_READ
            //IOobject::AUTO_WRITE
        ),
        mesh
     );
    const fvPatchField<scalar>& HCA = patch().lookupPatchField<volScalarField, scalar>(HCAName_);
    //scalarField&  HCA_ = HCA.ref();
    
    scalarField  ap = *this;
    tmp<scalarField> thetaCp(new scalarField(patch().size(),theta0_));
    //scalarField&  thetaC= thetaCp();
    scalarField&  thetaC = thetaCp.ref();    
    //scalarField thetaCp(patch().size(), 45.0);
    //scalarField thetaCp(size(), 45.0);
    //scalarField& thetaC = thetaCp;
   
    forAll(ap, faceI) //Loop trhough all faces on the alpha patch
    {
      thetaC[faceI] = HCA[faceI];    
    }

    return thetaCp; 
}


void Foam::heterogeneousAlphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("HCA", "HCA", HCAName_);
    os.writeEntry("theta0", theta0_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        heterogeneousAlphaContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
