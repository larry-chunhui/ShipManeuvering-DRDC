/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "twoDOscillatingDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    axis_(Zero),
    origin_(Zero),
    angle0_(0.0),
    amplitude_(0.0),
    omega_(0.0),
    p0_(p.localPoints()),
 //chliu   
    amplitudeH_(vector::zero),
    omegaH_(0.0)  
{}


twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    axis_(dict.lookup("axis")),
    origin_(dict.lookup("origin")),
    angle0_(readScalar(dict.lookup("angle0"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    omega_(readScalar(dict.lookup("omega"))),
//chliu
    amplitudeH_(dict.lookup("amplitudeH")),
    omegaH_(readScalar(dict.lookup("omegaH")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const twoDOscillatingDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    angle0_(ptf.angle0_),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_, mapper),
//chliu
    amplitudeH_(ptf.amplitudeH_),
    omegaH_(ptf.omegaH_)


{}


twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const twoDOscillatingDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    angle0_(ptf.angle0_),
    amplitude_(ptf.amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_),
//chliu
    amplitudeH_(ptf.amplitudeH_),
    omegaH_(ptf.omegaH_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void twoDOscillatingDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void twoDOscillatingDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const twoDOscillatingDisplacementPointPatchVectorField& aODptf =
        refCast<const twoDOscillatingDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void twoDOscillatingDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
//chunhui change from sin to cosine in 2018.10.29
    scalar angle = angle0_ + amplitude_*sin(omega_*t.value());
//    scalar angle = angle0_ + amplitude_*cos(omega_*t.value());

    vector axisHat = axis_/mag(axis_);
    vectorField p0Rel(p0_ - origin_);

if(t.value()<2.0*pi/omegaH_/4.0)
{    vectorField::operator=
    (
        p0Rel*(cos(angle) - 1)
      + (axisHat ^ p0Rel*sin(angle))
      + (axisHat & p0Rel)*(1 - cos(angle))*axisHat

    );
}
else
{
    vectorField::operator=
    (
//        p0Rel*(cos(angle) - 1)
//      + (axisHat ^ p0Rel*sin(angle))
//      + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
        p0Rel*(cos(angle) - 1)
      + (axisHat ^ p0Rel*sin(angle))
      + (-amplitudeH_*cos(omegaH_*t.value()))  //add chliu
      + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
    );
}
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void twoDOscillatingDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("axis", axis_);
    os.writeEntry("origin", origin_);
    os.writeEntry("angle0", angle0_);
    os.writeEntry("amplitude", amplitude_);
    os.writeEntry("omega", omega_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
//chliu
    os.writeKeyword("amplitudeH")
        << amplitudeH_ << token::END_STATEMENT << nl;
    os.writeKeyword("omegaH")
        << omegaH_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    twoDOscillatingDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
