/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "dissolutionPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(dissolutionPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        dissolutionPowerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::dissolutionPowerLaw::calcNu() const
{
    return (nu0_ - nuInf_)/(scalar(1) + pow(m_*strainRate(), n_)) + nuInf_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::dissolutionPowerLaw::dissolutionPowerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    dissolutionPowerLawCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    nu0_("nu0", dimViscosity, dissolutionPowerLawCoeffs_),
    nuInf_("nuInf", dimViscosity, dissolutionPowerLawCoeffs_),
    m_("m", dimTime, dissolutionPowerLawCoeffs_),
    
    nuSpe_("nuSpe", dimViscosity, dissolutionPowerLawCoeffs_),
    n_("n", dimless, dissolutionPowerLawCoeffs_),
    
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::dissolutionPowerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    dissolutionPowerLawCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    dissolutionPowerLawCoeffs_.readEntry("nu0", nu0_);
    dissolutionPowerLawCoeffs_.readEntry("nuInf", nuInf_);
    dissolutionPowerLawCoeffs_.readEntry("m", m_);
    
    dissolutionPowerLawCoeffs_.readEntry("nuSpe", nuSpe_);
    dissolutionPowerLawCoeffs_.readEntry("n", n_);

    return true;
}


// ************************************************************************* //
