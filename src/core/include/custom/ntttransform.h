#ifndef CUSTOM_NTTTRANSFORM_H
#define CUSTOM_NTTTRANSFORM_H

#include "math/hal/intnat/ubintnat.h"
#include "math/hal/intnat/mubintvecnat.h"

using namespace intnat;

struct FfiVec {
    FfiVec(unsigned long modulus,
        unsigned long* data,
        size_t size): m_modulus(modulus), m_data(data), m_size(size) {}
    unsigned long m_modulus;
    unsigned long* m_data;
    unsigned long m_size;
};

extern "C" {
    void rust_ntt_transform(
        const FfiVec* rootOfUnityTable,
        const FfiVec* preconRootOfUnityTable,
        FfiVec* element
    );
}

template <typename VecType>
void customNttWrapper(const VecType& rootOfUnityTable,
    const VecType& preconRootOfUnityTable,
    VecType* element) {
    FfiVec rootOfUnityTableWrap(static_cast<unsigned long>(rootOfUnityTable.GetModulus()),
     (unsigned long*)(&rootOfUnityTable[0]), rootOfUnityTable.GetLength());
    FfiVec preconRootOfUnityTableWrap(static_cast<unsigned long>(preconRootOfUnityTable.GetModulus()),
     (unsigned long*)(&preconRootOfUnityTable[0]), preconRootOfUnityTable.GetLength());
    FfiVec elementWrap(static_cast<unsigned long>(element->GetModulus()),
     (unsigned long*)(&(*element)[0]), element->GetLength());
    
    rust_ntt_transform(&rootOfUnityTableWrap, &preconRootOfUnityTableWrap, &elementWrap);
}

#endif
