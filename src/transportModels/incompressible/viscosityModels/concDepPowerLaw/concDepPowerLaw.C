/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd
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

\*---------------------------------------------------------------------------*/

#include "concDepPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(concDepPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        concDepPowerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::concDepPowerLaw::calcNu() const
{
    /*tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                U_.time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE  //AUTO_WRITE  NO_WRITE
            ),
        ),
        U_.mesh(),
        nuMax_
    );
    volScalarField& nu = tn.ref();*/
    
    // Check if the field exists before trying to access it
    if (U_.mesh().foundObject<volScalarField>(solutionSpeciesName_))
    {
        const volScalarField& Y_ = U_.mesh().lookupObject<volScalarField>(solutionSpeciesName_);
    
        return max
        (
            nuMin_,
            min
            (
                nuMax_,
                pow((pow(nuMax_,-0.25) + (Y_/YInitial_)* (pow(nuSpecies_,-0.25) - pow(nuMax_,-0.25))),-4)
            )
        );
    }
    else
    {
        tmp<volScalarField> tnu
	(
	    new volScalarField
	    (
		IOobject
		(
		    "nu",
		    U_.time().timeName(),
		    U_.mesh(),
		    IOobject::NO_READ,
		    IOobject::AUTO_WRITE  //AUTO_WRITE  NO_WRITE
		),
		U_.mesh(),
		nuMax_
	    )
	);
        return tnu;
    }
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::concDepPowerLaw::concDepPowerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    concDepPowerLawCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    k_("k", dimViscosity, concDepPowerLawCoeffs_),
    n_("n", dimless, concDepPowerLawCoeffs_),
    nuMin_("nuMin", dimViscosity, concDepPowerLawCoeffs_),
    nuMax_("nuMax", dimViscosity, concDepPowerLawCoeffs_),
    solutionSpeciesName_(concDepPowerLawCoeffs_.get<word>("solutionSpecies")),
    nuSpecies_("nuSpecies", dimViscosity, concDepPowerLawCoeffs_),
    YInitial_("C0",dimensionSet(0,-3,0,0,1,0,0),concDepPowerLawCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE  //AUTO_WRITE  NO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::concDepPowerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    concDepPowerLawCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    concDepPowerLawCoeffs_.readEntry("k", k_);
    concDepPowerLawCoeffs_.readEntry("n", n_);
    concDepPowerLawCoeffs_.readEntry("nuMin", nuMin_);
    concDepPowerLawCoeffs_.readEntry("nuMax", nuMax_);
    concDepPowerLawCoeffs_.readEntry("solutionSpecies", solutionSpeciesName_);
    concDepPowerLawCoeffs_.readEntry("nuSpecies", nuSpecies_);
    concDepPowerLawCoeffs_.readEntry("C0", YInitial_);
    

    return true;
}


// ************************************************************************* //
