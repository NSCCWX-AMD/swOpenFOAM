#include "Field.H"
#include "tensor.H"
#include "swField_host.h"

#if(SWTIMER)
//#include "Timers.hpp"
#endif


namespace Foam
{

//vector dot tensor

template< >
void dot <vector, tensor>
(Field<typename innerProduct<vector, tensor>::type>& res,
        const UList<vector>& f1, const UList<tensor>& f2)
{
    typedef typename innerProduct<vector, tensor>::type productType;

    checkFields(res, f1, f2, "f1 " "=" " f2 " "*" " f3");
    productType * __restrict__ f1P = (res).begin();
    const vector * __restrict__ f2P = (f1).begin();
    const tensor * __restrict__ f3P = (f2).begin();

    label i = (res).size();

#if 1
{
    printf("swField vector dot tensor call\n");
    printf("i=%ld\n", i);
    std::exit(0);
}
#endif

    swVectorDotTensor data;

    if (i > 2000) {
        data.size = i;
        data.vector_dim = 3;
        data.tensor_dim = 9;
        data.f1P = (scalar*) f1P;
        data.f2P = (scalar*) f2P;
        data.f3P = (scalar*) f3P;

        call_swVectorDotTensor_slave(&data);

    } else {
        while (i--) {
            (*f1P++) = (*f2P++) & (*f3P++);
        }
    }
}

//vector dot vector

template< >
void dot(Field<typename innerProduct<vector, vector>::type>& res, const UList<vector>& f1, const UList<vector>& f2)
{
    typedef typename innerProduct<vector, vector>::type productType;
    checkFields(res, f1, f2, "f1 " "=" " f2 " "&" " f3");
    productType * __restrict__ f1P = (res).begin();
    const vector * __restrict__ f2P = (f1).begin();
    const vector * __restrict__ f3P = (f2).begin();
    label i = (res).size();
    swVectorDotTensor data;

#if 1
{
    printf("swField vector dot vector call\n");
    printf("i=%ld\n", i);
    std::exit(0);
}
#endif

    if (i > 2000) {
        data.size = i;
        data.vector_dim = 3;

        data.f1P = (scalar*) f1P;
        data.f2P = (scalar*) f2P;
        data.f3P = (scalar*) f3P;

        call_swVectorDotVector_slave(&data);

    } else {
        while (i--) {
            (*f1P++) = (*f2P++) & (*f3P++);
        }
    }
}


////////////////////////////
template< >
void multiply(Field<scalar>& res, const scalar& s1, const UList<scalar>& f2)
{
    label i = (res).size();

    checkFields(res, f2, "f1 " "=" " s " "*" " f2");
    scalar * __restrict__ f1P = (res).begin();
    const scalar * __restrict__ f2P = (f2).begin();

    i = (res).size();
    if (i > 2000) {

        swVectorDotTensor data;
        data.size = i;
        data.f1P = (scalar*) f1P;
        data.f2P = (scalar*) f2P;
        data.s1 = s1;

        call_swVectorMulScalar_slave(&data);

    } else {
        while (i--) {
            (*f1P++) = (s1) * (*f2P++);
        }
    }
}

}//end namespace
