#ifndef SZ3_IMPL_SZDISPATCHER_HPP
#define SZ3_IMPL_SZDISPATCHER_HPP

#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/api/impl/SZInterp.hpp"
#include "QoZ/api/impl/SZLorenzoReg.hpp"
#include <cmath>


template<class T, QoZ::uint N>
char *SZ_compress_dispatcher(QoZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    QoZ::calAbsErrorBound(conf, data);

    char *cmpData;
    if (conf.cmprAlgo == QoZ::ALGO_LORENZO_REG) {
        cmpData = (char *) SZ_compress_LorenzoReg<T, N>(conf, data, outSize);
    } else if (conf.cmprAlgo == QoZ::ALGO_INTERP) {
        cmpData = (char *) SZ_compress_Interp<T, N>(conf, data, outSize);
    } else if (conf.cmprAlgo == QoZ::ALGO_INTERP_LORENZO) {
        cmpData = (char *) SZ_compress_Interp_lorenzo<T, N>(conf, data, outSize);
    }
    /*
    else if (conf.cmprAlgo == QoZ::ALGO_NEWINTERP) {
        cmpData = (char *) SZ_compress_NewInterp<T, N>(conf, data, outSize);
    }
    */
    else if (conf.cmprAlgo == QoZ::ALGO_INTERP_BLOCKED) {
        //std::cout<<"111"<<std::endl;
        cmpData = (char *) SZ_compress_Interp_blocked<T, N>(conf, data, outSize);
    }
    return cmpData;
}


template<class T, QoZ::uint N>
void SZ_decompress_dispatcher(QoZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    if (conf.cmprAlgo == QoZ::ALGO_LORENZO_REG) {
        SZ_decompress_LorenzoReg<T, N>(conf, cmpData, cmpSize, decData);
    } else if (conf.cmprAlgo == QoZ::ALGO_INTERP) {
        SZ_decompress_Interp<T, N>(conf, cmpData, cmpSize, decData);
    } else {
        printf("SZ_decompress_dispatcher, Method not supported\n");
        exit(0);
    }
}

#endif