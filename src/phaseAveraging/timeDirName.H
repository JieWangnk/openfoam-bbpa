/*---------------------------------------------------------------------------*\
Function
    Foam::timeDirName

Description
    The name of the time directory that holds the instant t, found by VALUE.

    A time name is a decimal rendering of a binary value. The spelling the
    solver wrote (at its precision, from its accumulated t) and the spelling
    Time::timeName(t) returns for the label computed here can differ in the
    last digits: on a production restart the bin products sat in
    "134.041276799808" while the restore asked for "134.04127679981", and the
    guard that refuses an unmatched checkpoint reset 44 of 48 bins. Matching
    the directory by its numeric value within tol removes the dependence on
    spelling.

    The exact spelling is tried first (the common case, one stat call). With
    no directory within tol on disk the exact spelling is returned, so a
    writer about to create a new directory behaves exactly as before. The
    "constant" directory is never a candidate.

    tol must stay below half the solver step, so that at most one solver
    write directory can lie within it; callers pass
    min(1e-6*binDeltaT, 0.25*deltaT).

\*---------------------------------------------------------------------------*/

#ifndef timeDirName_H
#define timeDirName_H

#include "Time.H"
#include "instantList.H"
#include "OSspecific.H"

namespace Foam
{

inline word timeDirName(const Time& time, const scalar t, const scalar tol)
{
    const word exact = Time::timeName(t);
    if (isDir(time.path()/exact)) return exact;

    const instantList times = time.times();
    label best = -1;
    scalar dBest = tol;
    forAll(times, i)
    {
        if (times[i].name() == time.constant()) continue;
        const scalar d = mag(times[i].value() - t);
        if (d <= dBest)
        {
            dBest = d;
            best = i;
        }
    }
    return best >= 0 ? times[best].name() : exact;
}


//- The lookup tolerance for a bin of width binDeltaT on a run stepping at dt
inline scalar timeDirTol(const scalar binDeltaT, const scalar dt)
{
    return min(1e-6*binDeltaT, 0.25*dt);
}

} // End namespace Foam

#endif
