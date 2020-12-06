#ifndef __Remnants_h__
#define __Remnants_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeGB.h"


class BaseStar;
class HeGB;

class Remnants: virtual public BaseStar, public HeGB {

public:

    Remnants(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeGB(p_BaseStar, false) {
        if (p_Initialise) Initialise();
    }

    Remnants& operator = (const BaseStar &p_BaseStar) { static_cast<BaseStar&>(*this) = p_BaseStar; return *this; }


    // member functions


protected:


    // member functions - alphabetically

};

#endif // __Remnants_h__