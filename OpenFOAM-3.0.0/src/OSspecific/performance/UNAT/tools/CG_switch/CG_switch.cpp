extern "C"
{
#include "athread_switch.h"
}

namespace Foam
{

// class that initialize and
// halt slave cores before and
// after main function
class CG_switch
{
    public:
        CG_switch()
        {
            CG_init();
        }
        ~CG_switch()
        {
            CG_halt();
        }
};

CG_switch cg_switch;

}//namespace Foam
