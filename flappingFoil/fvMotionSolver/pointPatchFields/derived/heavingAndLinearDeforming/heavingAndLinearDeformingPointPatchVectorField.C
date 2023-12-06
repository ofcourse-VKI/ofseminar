/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "heavingAndLinearDeformingPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heavingAndLinearDeformingPointPatchVectorField::
heavingAndLinearDeformingPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    ha_(Zero),
    frequency_(0.0),
    hp_(Zero),
    xoff_(0.05),
    a0_(Zero),
    a1_(Zero),
    length_(0.0)
{}


Foam::heavingAndLinearDeformingPointPatchVectorField::
heavingAndLinearDeformingPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    ha_(dict.lookup("ha")),
    frequency_(dict.getOrDefault<scalar>("frequency", Zero)),
    hp_(dict.lookup("hp")),
    xoff_(dict.get<scalar>("xoff")),
    a0_(dict.lookup("a0")),
    a1_(dict.lookup("a1")),
    length_(dict.get<scalar>("length"))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::heavingAndLinearDeformingPointPatchVectorField::
heavingAndLinearDeformingPointPatchVectorField
(
    const heavingAndLinearDeformingPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    ha_(ptf.ha_),
    frequency_(ptf.frequency_),
    hp_(ptf.hp_),
    xoff_(ptf.xoff_),
    a0_(ptf.a0_),
    a1_(ptf.a1_),
    length_(ptf.length_)
{}


Foam::heavingAndLinearDeformingPointPatchVectorField::
heavingAndLinearDeformingPointPatchVectorField
(
    const heavingAndLinearDeformingPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    ha_(ptf.ha_),
    frequency_(ptf.frequency_),
    hp_(ptf.hp_),
    xoff_(ptf.xoff_),
    a0_(ptf.a0_),
    a1_(ptf.a1_),
    length_(ptf.length_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::heavingAndLinearDeformingPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    //const scalarField points( offset_ & patch().localPoints());
    scalarField xCoord = patch().localPoints().component(vector::X);
    
    
    Field<vector>::operator=
    (
        (ha_*cos(2*M_PI*frequency_*t.value()) - ha_) +
	 hp_*sin(2*M_PI*frequency_*t.value())*sqr(xCoord-xoff_)
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::heavingAndLinearDeformingPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("ha", ha_);
    os.writeEntry("frequency", frequency_);
    os.writeEntry("xoff", xoff_);
    os.writeEntry("hp", hp_);
    os.writeEntry("a0", a0_);
    os.writeEntry("a1", a1_);
    os.writeEntry("length", length_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        heavingAndLinearDeformingPointPatchVectorField
    );
}

// ************************************************************************* //
