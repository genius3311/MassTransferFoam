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

Class
    Foam::heterogeneousAlphaContactAngleFvPatchScalarField

Description
    A constant alphaContactAngle scalar boundary condition
    (alphaContactAngleTwoPhaseFvPatchScalarField)

SourceFiles
    heterogeneousAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef heterogeneousAlphaContactAngleFvPatchScalarField_H
#define heterogeneousAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class constantAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class heterogeneousAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Name of temperature field, default = "T"
        word HCAName_;

        //- Equilibrium contact angle function
        scalar theta0_;
        


public:

    //- Runtime type information
    TypeName("heterogeneousAlphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        heterogeneousAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        heterogeneousAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  heterogeneousAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        heterogeneousAlphaContactAngleFvPatchScalarField
        (
            const heterogeneousAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        heterogeneousAlphaContactAngleFvPatchScalarField
        (
            const heterogeneousAlphaContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
       /* virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new heterogeneousAlphaContactAngleFvPatchScalarField(*this)
            );
        }*/

        //- Construct as copy setting internal field reference
        heterogeneousAlphaContactAngleFvPatchScalarField
        (
            const heterogeneousAlphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        /*virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new heterogeneousAlphaContactAngleFvPatchScalarField(*this, iF)
            );
        }*/
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
