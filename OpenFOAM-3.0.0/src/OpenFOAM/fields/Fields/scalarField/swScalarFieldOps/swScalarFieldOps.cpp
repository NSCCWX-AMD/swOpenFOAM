/*---------------------------------------------------------------------------*\
Copyright (C) 2011 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of CAELUS.

    CAELUS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CAELUS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CAELUS.  If not, see <http://www.gnu.org/licenses/>.

Description
    Operations' specialisation of Field\<T\> for scalar.

\*---------------------------------------------------------------------------*/

#include "swScalarFieldOps.hpp"
#include "vectorOps.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//template <>
void multiply
(
    Field<scalar>& res,
    const UList<scalar>& f1,
    const UList<scalar>& f2
)
{
    checkFields(res, f1, f2, "f1 " "=" " f2 " "*" " f3");
    scalar* __restrict__ f1P = (res).begin();
    const scalar* __restrict__ f2P = (f1).begin();
    const scalar* __restrict__ f3P = (f2).begin();
    label i = (res).size();

#if 1
{
    printf("swScalarFieldOps multiply ncall\n");
    printf("accUsingSize=%d\n", accUsingSize);
    std::exit(0);
}
#endif

    // test the slave function
    if(i > accUsingSize)
    {
#if 1
{
    static int ncall = 0;
    ncall++;
    if(ncall%1==0)
    {
        printf("swScalarFieldOps multiply ncall=%d\n", ncall);
    }
}
#endif
        MVM_Arrays paras;
        init_MVM_Arrays(&paras, i);
        paras.A1Ptr = (scalar*) f1P;
        paras.A2Ptr = (scalar*) f2P;
        paras.A3Ptr = (scalar*) f3P;
        vectorOps_host(&paras, &slave_userFunc_aEbMuc);
    }
    else
    {
        while (i--)
        {
            (*f1P++) = (*f2P++) * (*f3P++);
        }
    }

}

void divide
(
    Field<scalar>& res,
    const UList<scalar>& f1,
    const UList<scalar>& f2
)
{
    checkFields(res, f1, f2, "f1 " "=" " f2 " "/" " f3");
    scalar* __restrict__ f1P = (res).begin();
    const scalar* __restrict__ f2P = (f1).begin();
    const scalar* __restrict__ f3P = (f2).begin();
    label i = (res).size();

#if 1
{
    printf("swScalarFieldOps divide ncall\n");
    printf("accUsingSize=%d\n", accUsingSize);
    std::exit(0);
}
#endif

    // test the slave function
    if(i > accUsingSize)
    {
        MVM_Arrays paras;
        init_MVM_Arrays(&paras, i);
        paras.A1Ptr = (scalar*) f1P;
        paras.A2Ptr = (scalar*) f2P;
        paras.A3Ptr = (scalar*) f3P;
        vectorOps_host(&paras, &slave_userFunc_aEbDc);
    }
    else
    {
        while (i--)
        {
            (*f1P++) = (*f2P++) / (*f3P++);
        }
    }

}
} // End namespace CML
