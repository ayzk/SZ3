#ifndef SZ3_SZ_LORENZO_REG_HPP
#define SZ3_SZ_LORENZO_REG_HPP

#include "QoZ/compressor/SZGeneralCompressor.hpp"
#include "QoZ/frontend/SZFastFrontend.hpp"
#include "QoZ/frontend/SZGeneralFrontend.hpp"
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/predictor/ComposedPredictor.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/predictor/RegressionPredictor.hpp"
#include "QoZ/predictor/PolyRegressionPredictor.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
#include "QoZ/utils/Iterator.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Extraction.hpp"
#include "QoZ/utils/QuantOptimizatioin.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/def.hpp"
#include <cmath>
#include <memory>


template<class T, QoZ::uint N, class Quantizer, class Encoder, class Lossless>
std::shared_ptr<QoZ::concepts::CompressorInterface<T>>
make_lorenzo_regression_compressor(const QoZ::Config &conf, Quantizer quantizer, Encoder encoder, Lossless lossless) {
    std::vector<std::shared_ptr<QoZ::concepts::PredictorInterface<T, N>>> predictors;

    int methodCnt = (conf.lorenzo + conf.lorenzo2 + conf.regression + conf.regression2);
    int use_single_predictor = (methodCnt == 1);
    if (methodCnt == 0) {
        printf("All lorenzo and regression methods are disabled.\n");
        exit(0);
    }
    if (conf.lorenzo) {
        std::vector<double> coeffs;
        if(conf.useCoeff)
            coeffs=conf.lorenzo1_coeffs;
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound,coeffs), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound,coeffs));
        }
    }
    if (conf.lorenzo2) {
        std::vector<double> coeffs;
        if(conf.useCoeff)
            coeffs=conf.lorenzo2_coeffs;
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound,coeffs), quantizer),
                    encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound,coeffs));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::RegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::RegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }

    if (conf.regression2) {
        if (use_single_predictor) {
            return QoZ::make_sz_general_compressor<T, N>(
                    QoZ::make_sz_general_frontend<T, N>(conf, QoZ::PolyRegressionPredictor<T, N>(conf.blockSize, conf.absErrorBound),
                                                       quantizer), encoder, lossless);
        } else {
            predictors.push_back(std::make_shared<QoZ::PolyRegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }
    return QoZ::make_sz_general_compressor<T, N>(
            QoZ::make_sz_general_frontend<T, N>(conf, QoZ::ComposedPredictor<T, N>(predictors),
                                               quantizer), encoder, lossless);
}


template<class T, QoZ::uint N>
char *SZ_compress_LorenzoReg(QoZ::Config &conf, T *data, size_t &outSize) {

    assert(N == conf.N);
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);
    //QoZ::calAbsErrorBound(conf, data);

    char *cmpData;
    auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
    if (N == 3 and !conf.regression2 and !conf.useCoeff) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer), QoZ::HuffmanEncoder<int>(),
                                                       QoZ::Lossless_zstd());
        cmpData = (char *) sz->compress(conf, data, outSize);
    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        //std::cout<<"lor1"<<std::endl;
        cmpData = (char *) sz->compress(conf, data, outSize);
    }
    return cmpData;
}


template<class T, QoZ::uint N>
void SZ_decompress_LorenzoReg(const QoZ::Config &conf, char *cmpData, size_t cmpSize, T *decData) {
    assert(conf.cmprAlgo == QoZ::ALGO_LORENZO_REG);

    QoZ::uchar const *cmpDataPos = (QoZ::uchar *) cmpData;
    QoZ::LinearQuantizer<T> quantizer;
    if (N == 3 and !conf.regression2 and !conf.useCoeff) {
        // use fast version for 3D
        auto sz = QoZ::make_sz_general_compressor<T, N>(QoZ::make_sz_fast_frontend<T, N>(conf, quantizer),
                                                       QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);
        return;

    } else {
        auto sz = make_lorenzo_regression_compressor<T, N>(conf, quantizer, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        sz->decompress(cmpDataPos, cmpSize, decData);
        return;
    }

}

#endif