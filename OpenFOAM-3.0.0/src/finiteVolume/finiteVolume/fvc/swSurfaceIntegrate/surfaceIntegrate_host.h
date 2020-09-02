#include<stdio.h>
#ifndef SURFACEINTEGRATE_HOST_H
#define SURFACEINTEGRATE_HOST_H

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{

 void surfaceIntegrate_host(const int*,const int*, const double*,double*,int,int,int,int,int);

}

}
#endif
