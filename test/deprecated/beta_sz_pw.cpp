#include "QoZ/frontend/SZGeneralFrontend.hpp"
#include "QoZ/predictor/Predictor.hpp"
#include "QoZ/predictor/LorenzoPredictor.hpp"
#include "QoZ/predictor/RegressionPredictor.hpp"
#include "QoZ/predictor/ComposedPredictor.hpp"
#include "QoZ/quantizer/IntegerQuantizer.hpp"
#include "QoZ/utils/FileUtil.hpp"
#include "QoZ/utils/Config.hpp"
#include "QoZ/def.hpp"
#include "QoZ/compressor/SZGeneralCompressor.hpp"
#include "QoZ/encoder/HuffmanEncoder.hpp"
#include "QoZ/encoder/ArithmeticEncoder.hpp"
#include "QoZ/encoder/BypassEncoder.hpp"
#include "QoZ/lossless/Lossless_zstd.hpp"
#include "QoZ/lossless/Lossless_bypass.hpp"
#include "QoZ/utils/Statistic.hpp"
#include "QoZ/utils/Timer.hpp"
#include "QoZ/version.hpp"
#include <sstream>
#include <random>
#include <cstdio>
#include <iostream>
#include <memory>

std::string src_file_name;
float relative_error_bound = 0;
float compression_time = 0;


template<typename T, uint N, class Frontend, class Encoder, class Lossless>
float SZ_compress(std::unique_ptr<T[]> const &data,
                  QoZ::Config &conf,
                  Frontend frontend, Encoder encoder, Lossless lossless) {

    std::cout << "****************** Options ********************" << std::endl;
    std::cout << "dimension = " << N
              << ", error bound = " << conf.absErrorBound
              << ", blockSize = " << conf.blockSize
              << ", stride = " << conf.stride
              << ", quan_state_num = " << conf.quantbinCnt
              << std::endl
              << "lorenzo = " << conf.lorenzo
              << ", 2ndlorenzo = " << conf.lorenzo2
              << ", regression = " << conf.regression
              << ", encoder = " << conf.encoder
              << ", lossless = " << conf.lossless
              << std::endl;


    std::vector<T> data_ = std::vector<T>(data.get(), data.get() + conf.num);
    // preprocess
    std::vector<unsigned char> signs(conf.num, 0);
    float max_abs_log_data = 0;
    float min_log_data = std::numeric_limits<float>::max();
    bool positive = true;
    for (size_t i = 0; i < conf.num; i++) {
        if (data[i] == 0) {
            signs[i] = 2;
        }
        if (data[i] < 0) {
            signs[i] = 1;
            data[i] = -data[i];
            positive = false;
        }
        if (data[i] > 0) {
            data[i] = log2(data[i]);
            if (data[i] > max_abs_log_data) max_abs_log_data = data[i];
            if (data[i] < min_log_data) min_log_data = data[i];
        }

    }

    if (fabs(min_log_data) > max_abs_log_data) max_abs_log_data = fabs(min_log_data);
    double realPrecision = log2(1.0 + conf.absErrorBound) - max_abs_log_data * 1.2e-7;
    conf.absErrorBound = realPrecision;
    for (size_t i = 0; i < conf.num; i++) {
        if (signs[i] == 2) {
            data[i] = min_log_data - 2.0001 * realPrecision;
            signs[i] = 0;
        }
    }

    double minLogValue = min_log_data - 1.0001 * realPrecision;
    std::cout << "log eb=" << conf.absErrorBound
              << ", minlogvalue=" << minLogValue << std::endl;
    if (!positive) {
        // compress signs
//        unsigned long signSize = sz_lossless_compress(ZSTD_COMPRESSOR, 3, signs, dataLength, &comp_signs);
    }

//    std::vector<T> data1 = std::vector<T>(data.get(), data.get() + conf.num);
    auto sz = QoZ::make_sz_general_compressor<T, N>(frontend, encoder, lossless);

    QoZ::Timer timer;
    timer.start();
    std::cout << "****************** Compression ******************" << std::endl;


    size_t compressed_size = 0;
    std::unique_ptr<QoZ::uchar[]> compressed;
    compressed.reset(sz->compress(conf, data.get(), compressed_size));

    compression_time = timer.stop("Compression");

    auto ratio = conf.num * sizeof(T) * 1.0 / compressed_size;
    std::cout << "Compression Ratio = " << ratio << std::endl;
    std::cout << "Compressed size = " << compressed_size << std::endl;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 10000);
    std::stringstream ss;
    ss << src_file_name.substr(src_file_name.rfind('/') + 1)
       << "." << relative_error_bound << "." << dis(gen) << ".sz3";
    auto compressed_file_name = ss.str();
    QoZ::writefile(compressed_file_name.c_str(), compressed.get(), compressed_size);
    std::cout << "Compressed file = " << compressed_file_name << std::endl;

    std::cout << "****************** Decompression ****************" << std::endl;
    compressed = QoZ::readfile<QoZ::uchar>(compressed_file_name.c_str(), compressed_size);

    timer.start();
    T *dec_data = sz->decompress(compressed.get(), compressed_size, conf.num);
    timer.stop("Decompression");
//    QoZ::verify<T>(data1.data(), dec_data, conf.num);

//    float threshold = minLogValue;
    if (!positive) {
//        sz_lossless_decompress(ZSTD_COMPRESSOR, tdps->pwrErrBoundBytes, tdps->pwrErrBoundBytes_size, &signs, dataSeriesLength);
        for (size_t i = 0; i < conf.num; i++) {
            if (dec_data[i] < minLogValue) dec_data[i] = 0;
            else dec_data[i] = exp2(dec_data[i]);
            if (signs[i]) dec_data[i] = -dec_data[i];
        }
    } else {
        for (size_t i = 0; i < conf.num; i++) {
            if (dec_data[i] < minLogValue) dec_data[i] = 0;
            else dec_data[i] = exp2(dec_data[i]);
        }
    }

    QoZ::verify<T>(data_.data(), dec_data, conf.num);

    remove(compressed_file_name.c_str());
    auto decompressed_file_name = compressed_file_name + ".out";
    QoZ::writefile(decompressed_file_name.c_str(), dec_data, conf.num);
    std::cout << "Decompressed file = " << decompressed_file_name << std::endl;

    delete[] dec_data;
    return ratio;
}

template<typename T, uint N, class Frontend>
float SZ_compress_build_backend(std::unique_ptr<T[]> const &data,
                                QoZ::Config &conf,
                                Frontend frontend) {
    if (conf.lossless == 1) {
        if (conf.encoder == 1) {
            return SZ_compress<T, N>(data, conf, frontend, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_zstd());
        } else if (conf.encoder == 2) {
            return SZ_compress<T, N>(data, conf, frontend, QoZ::ArithmeticEncoder<int>(true), QoZ::Lossless_zstd());
        } else {
            return SZ_compress<T, N>(data, conf, frontend, QoZ::BypassEncoder<int>(), QoZ::Lossless_zstd());
        }
    } else {
        if (conf.encoder == 1) {
            return SZ_compress<T, N>(data, conf, frontend, QoZ::HuffmanEncoder<int>(), QoZ::Lossless_bypass());
        } else if (conf.encoder == 2) {
            return SZ_compress<T, N>(data, conf, frontend, QoZ::ArithmeticEncoder<int>(true), QoZ::Lossless_bypass());
        } else {
            return SZ_compress<T, N>(data, conf, frontend, QoZ::BypassEncoder<int>(), QoZ::Lossless_bypass());
        }
    }
}

template<typename T, uint N>
float SZ_compress_build_frontend(std::unique_ptr<T[]> const &data, QoZ::Config &conf) {
    auto quantizer = QoZ::LinearQuantizer<T>(conf.absErrorBound, conf.quantbinCnt / 2);
    std::vector<std::shared_ptr<QoZ::concepts::PredictorInterface<T, N>>> predictors;

    int use_single_predictor =
            (conf.lorenzo + conf.lorenzo2 + conf.regression) == 1;
    if (conf.lorenzo) {
        if (use_single_predictor) {
            return SZ_compress_build_backend<T, N>(data, conf,
                                                   QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 1>(conf.absErrorBound),
                                                                              quantizer));
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 1>>(conf.absErrorBound));
        }
    }
    if (conf.lorenzo2) {
        if (use_single_predictor) {
            return SZ_compress_build_backend<T, N>(data, conf,
                                                   QoZ::make_sz_general_frontend<T, N>(conf, QoZ::LorenzoPredictor<T, N, 2>(conf.absErrorBound),
                                                                              quantizer));
        } else {
            predictors.push_back(std::make_shared<QoZ::LorenzoPredictor<T, N, 2>>(conf.absErrorBound));
        }
    }
    if (conf.regression) {
        if (use_single_predictor) {
            return SZ_compress_build_backend<T, N>(data, conf,
                                                   QoZ::make_sz_general_frontend<T, N>(conf, QoZ::RegressionPredictor<T, N>(conf.blockSize,
                                                                                                                  conf.absErrorBound),
                                                                              quantizer));
        } else {
            predictors.push_back(std::make_shared<QoZ::RegressionPredictor<T, N>>(conf.blockSize, conf.absErrorBound));
        }
    }

    return SZ_compress_build_backend<T, N>(data, conf,
                                           QoZ::make_sz_general_frontend<T, N>(conf, QoZ::ComposedPredictor<T, N>(predictors), quantizer));
}


template<class T, uint N>
float SZ_compress_parse_args(int argc, char **argv, int argp, std::unique_ptr<T[]> &data, float eb,
                             std::array<size_t, N> dims) {
    QoZ::Config conf;
    conf.setDims(dims.begin(), dims.end());
    conf.absErrorBound = eb;
    if (argp < argc) {
        int block_size = atoi(argv[argp++]);
        conf.blockSize = block_size;
        conf.stride = block_size;
    }

    int lorenzo_op = 1;
    if (argp < argc) {
        lorenzo_op = atoi(argv[argp++]);
        conf.lorenzo = lorenzo_op == 1 || lorenzo_op == 3;
        conf.lorenzo2 = lorenzo_op == 2 || lorenzo_op == 3;
    }

    int regression_op = 1;
    if (argp < argc) {
        regression_op = atoi(argv[argp++]);
        conf.regression = regression_op == 1;
    }

    if (argp < argc) {
        conf.encoder = atoi(argv[argp++]);
        if (conf.encoder == 0) {
            conf.quantbinCnt = 250;
        } else if (conf.encoder == 2) {
            conf.quantbinCnt = 1024;
        }
    }

    if (argp < argc) {
        conf.lossless = atoi(argv[argp++]);
    }

    if (argp < argc) {
        conf.quantbinCnt = atoi(argv[argp++]);
    }

    auto ratio = SZ_compress_build_frontend<T, N>(data, conf);
    return ratio;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "SZ v" << PROJECT_VER << std::endl;
        std::cout << "usage: " << argv[0] <<
                  " data_file -num_dim dim0 .. dimn relative_eb [blocksize lorenzo_op regression_op encoder lossless quantbinCnt]"
                  << std::endl;
        std::cout << "example: " << argv[0] <<
                  " qmcpack.dat -3 33120 69 69 1e-3 [6 1 1 1 1 32768]" << std::endl;
        return 0;
    }

    size_t num = 0;
    auto data = QoZ::readfile<float>(argv[1], num);
    src_file_name = argv[1];
    std::cout << "Read " << num << " elements\n";

    int dim = atoi(argv[2] + 1);
    assert(1 <= dim && dim <= 4);
    int argp = 3;
    std::vector<size_t> dims(dim);
    for (int i = 0; i < dim; i++) {
        dims[i] = atoi(argv[argp++]);
    }

//    char *eb_op = argv[argp++] + 1;
    float eb = 0;
//    if (*eb_op == 'a') {
//        eb = atof(argv[argp++]);
//    } else {
    eb = atof(argv[argp++]);
//    float max = data[0];
//    float min = data[0];
//    for (int i = 1; i < num; i++) {
//        if (max < data[i]) max = data[i];
//        if (min > data[i]) min = data[i];
//    }
//    eb = relative_error_bound * (max - min);
//    }

    if (dim == 1) {
        SZ_compress_parse_args<float, 1>(argc, argv, argp, data, eb, std::array<size_t, 1>{dims[0]});
    } else if (dim == 2) {
        SZ_compress_parse_args<float, 2>(argc, argv, argp, data, eb, std::array<size_t, 2>{dims[0], dims[1]});
    } else if (dim == 3) {
        SZ_compress_parse_args<float, 3>(argc, argv, argp, data, eb,
                                         std::array<size_t, 3>{dims[0], dims[1], dims[2]});
    } else if (dim == 4) {
        SZ_compress_parse_args<float, 4>(argc, argv, argp, data, eb,
                                         std::array<size_t, 4>{dims[0], dims[1], dims[2], dims[3]});
    }


    return 0;
}
