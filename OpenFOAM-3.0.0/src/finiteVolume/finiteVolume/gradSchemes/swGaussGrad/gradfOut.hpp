#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"
#include <fstream>

#ifndef gradfOut_HPP
#define gradfOut_HPP


namespace Foam
{

template<class Type, class GradType>
void gradfOut(const vectorField &Sf, const Field<Type>& issf, const Field<GradType>& igGrad)
{
    Info<<"gradfOut is not implemented!!"<<endl;
    ::exit(-1);
}


template<>
void gradfOut< vector , tensor>(
        const vectorField &Sf,
        const Field<vector>& issf,
        const Field<tensor>& igGrad
    )
{
#if 1
    std::ofstream faceOut("faceOut");
    std::ofstream cellOut("cellOut");
    for( label i = 0; i < Sf.size(); i++){
        faceOut <<Sf.cdata()[i].v_[0]<<'\t'
                <<Sf.cdata()[i].v_[1]<<'\t'
                <<Sf.cdata()[i].v_[2]<<"; ";
        faceOut <<issf.cdata()[i].v_[0]<<'\t'
                <<issf.cdata()[i].v_[1]<<'\t'
                <<issf.cdata()[i].v_[2]<<"; "
                <<"\n";
    }
    for( label i = 0; i < igGrad.size(); i++){
        cellOut <<igGrad.cdata()[i].v_[0]<<'\t'
                <<igGrad.cdata()[i].v_[1]<<'\t'
                <<igGrad.cdata()[i].v_[2]<<'\t'
                <<igGrad.cdata()[i].v_[3]<<'\t'
                <<igGrad.cdata()[i].v_[4]<<'\t'
                <<igGrad.cdata()[i].v_[5]<<'\t'
                <<igGrad.cdata()[i].v_[6]<<'\t'
                <<igGrad.cdata()[i].v_[7]<<'\t'
                <<igGrad.cdata()[i].v_[8]<<"; "
                <<"\n";
    }
    faceOut.close();
    cellOut.close();
    ::exit(0);
#endif
}

template<>
void gradfOut< scalar , vector >(
        const vectorField &Sf,
        const Field<scalar>& issf,
        const Field<vector>& igGrad
    )
{
#if 0
    std::ofstream faceOut("faceOut");
    std::ofstream cellOut("cellOut");
    for( label i = 0; i < Sf.size(); i++){
        faceOut <<Sf.cdata()[i].v_[0]<<'\t'
                <<Sf.cdata()[i].v_[1]<<'\t'
                <<Sf.cdata()[i].v_[2]<<"; ";
        faceOut <<issf.cdata()[i]<<"; "
                <<"\n";
    }
    for( label i = 0; i < igGrad.size(); i++){
        cellOut <<igGrad.cdata()[i].v_[0]<<'\t'
                <<igGrad.cdata()[i].v_[1]<<'\t'
                <<igGrad.cdata()[i].v_[2]<<"; "
                <<"\n";
    }
    faceOut.close();
    cellOut.close();
    ::exit(0);
#endif
}

}//namespace Foam

#endif
