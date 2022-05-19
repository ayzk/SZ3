//
// Created by Kai Zhao on 1/28/21.
//

#ifndef SZ3_BYTEUTIL_HPP
#define SZ3_BYTEUTIL_HPP

#include "def.hpp"
#include "MemoryUtil.hpp"
#include <cstring>

namespace SZ {

    typedef union lint16 {
        unsigned short usvalue;
        short svalue;
        unsigned char byte[2];
    } lint16;

    typedef union lint32 {
        int ivalue;
        unsigned int uivalue;
        unsigned char byte[4];
    } lint32;

    typedef union lint64 {
        long lvalue;
        unsigned long ulvalue;
        unsigned char byte[8];
    } lint64;

    typedef union ldouble {
        double value;
        unsigned long lvalue;
        unsigned char byte[8];
    } ldouble;

    typedef union lfloat {
        float value;
        unsigned int ivalue;
        unsigned char byte[4];
        uint16_t int16[2];
    } lfloat;

    inline void symTransform_4bytes(uchar data[4]) {
        unsigned char tmp = data[0];
        data[0] = data[3];
        data[3] = tmp;

        tmp = data[1];
        data[1] = data[2];
        data[2] = tmp;
    }

    inline unsigned short bytesToUInt16_bigEndian(const uchar *bytes) {
        int temp = 0;
        unsigned short res = 0;

        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        return res;
    }

    inline unsigned int bytesToUInt32_bigEndian(const uchar *bytes) {
        unsigned int temp = 0;
        unsigned int res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline unsigned long bytesToUInt64_bigEndian(const uchar *b) {
        unsigned long temp = 0;
        unsigned long res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }

    inline short bytesToInt16_bigEndian(unsigned char *bytes) {
        int temp = 0;
        short res = 0;

        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        return res;
    }

    inline int bytesToInt32_bigEndian(unsigned char *bytes) {
        int temp = 0;
        int res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline long bytesToInt64_bigEndian(const unsigned char *b) {
        long temp = 0;
        long res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }

    inline int bytesToInt_bigEndian(const unsigned char *bytes) {
        int temp = 0;
        int res = 0;

        res <<= 8;
        temp = bytes[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = bytes[3] & 0xff;
        res |= temp;

        return res;
    }

    inline long bytesToLong_bigEndian(const uchar *b) {
        long temp = 0;
        long res = 0;

        res <<= 8;
        temp = b[0] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[1] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[2] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[3] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[4] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[5] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[6] & 0xff;
        res |= temp;

        res <<= 8;
        temp = b[7] & 0xff;
        res |= temp;

        return res;
    }

    inline void longToBytes_bigEndian(unsigned char *b, unsigned long num) {
        b[0] = (unsigned char) (num >> 56);
        b[1] = (unsigned char) (num >> 48);
        b[2] = (unsigned char) (num >> 40);
        b[3] = (unsigned char) (num >> 32);
        b[4] = (unsigned char) (num >> 24);
        b[5] = (unsigned char) (num >> 16);
        b[6] = (unsigned char) (num >> 8);
        b[7] = (unsigned char) (num);
    }

    inline void int16ToBytes_bigEndian(unsigned char *b, uint16_t num) {
        b[0] = (unsigned char) (num >> 8);
        b[1] = (unsigned char) (num);
    }

    inline void int32ToBytes_bigEndian(unsigned char *b, uint32_t num) {
        b[0] = (unsigned char) (num >> 24);
        b[1] = (unsigned char) (num >> 16);
        b[2] = (unsigned char) (num >> 8);
        b[3] = (unsigned char) (num);
    }

    inline void intToBytes_bigEndian(unsigned char *b, unsigned int num) {
        b[0] = (unsigned char) (num >> 24);
        b[1] = (unsigned char) (num >> 16);
        b[2] = (unsigned char) (num >> 8);
        b[3] = (unsigned char) (num);
    }

    inline void int64ToBytes_bigEndian(unsigned char *b, uint64_t num) {
        b[0] = (unsigned char) (num >> 56);
        b[1] = (unsigned char) (num >> 48);
        b[2] = (unsigned char) (num >> 40);
        b[3] = (unsigned char) (num >> 32);
        b[4] = (unsigned char) (num >> 24);
        b[5] = (unsigned char) (num >> 16);
        b[6] = (unsigned char) (num >> 8);
        b[7] = (unsigned char) (num);
    }

    int bytesToInt(const unsigned char *bytes) {
        lfloat buf;
        memcpy(buf.byte, bytes, 4);
        return buf.ivalue;
    }

    std::string floatToBinary(float f) {
        lfloat u;
        u.value = f;
        std::string str(32, '0');
        for (int i = 0; i < 32; i++) {
            str[31 - i] = (u.ivalue % 2) ? '1' : '0';
            u.ivalue >>= 1;
        }
        return str;
    }

    template<class T>
    void truncateArray(T data, size_t n, int byteLen, uchar *&binary) {
        lfloat bytes;
        int b;
        for (size_t i = 0; i < n; i++) {
            bytes.value = data[i];
            for (b = 4 - byteLen; b < 4; b++) {
                *binary++ = bytes.byte[b];
            }
//            std::cout << std::bitset<32>(data[i]) << " " << std::bitset<16>(*binary) << '\n';
        }
    }

    template<class T>
    void truncateArrayRecover(uchar *binary, size_t n, int byteLen, T *data) {
        lfloat bytes;
        bytes.ivalue = 0;
        int b;
        for (size_t i = 0; i < n; i++) {
            for (b = 4 - byteLen; b < 4; b++) {
                bytes.byte[b] = *binary++;
            }
            data[i] = bytes.value;
        }
    }

    std::vector<uchar> LeadingBitsEncode(float pre, float data) {
        lfloat lfBuf_pre;
        lfloat lfBuf_cur;

        lfBuf_pre.value = pre;
        lfBuf_cur.value = data;
        lfBuf_pre.ivalue = lfBuf_cur.ivalue ^ lfBuf_pre.ivalue;

        std::vector<uchar> bytes;
        int n = 0;
        if (lfBuf_pre.ivalue == 0) {
            n = 0;
        } else if (lfBuf_pre.ivalue >> 8 == 0) {
            n = 1;
        } else if (lfBuf_pre.ivalue >> 16 == 0) {
            n = 2;
        } else if (lfBuf_pre.ivalue >> 24 == 0) {
            n = 3;
        } else {
            n = 4;
        }

        for (int i = 0; i < n; i++) {
            bytes.push_back(lfBuf_cur.byte[i]);
        }
        return bytes;
    }

    float LeadingBitsDecode(float pre, std::vector<uchar> bytes) {
        lfloat lfBuf_pre;
        lfloat lfBuf_cur;

        lfBuf_pre.value = pre;
        lfBuf_cur = lfBuf_pre;

        for (int i = 0; i < bytes.size(); i++) {
            lfBuf_cur.byte[i] = bytes[i];
        }
        return lfBuf_cur.value;
    }

    inline void encode_int_2bits(const std::vector<int> &data, uchar *&c) {

        size_t intLen = data.size();
        size_t byteLen = intLen * 2 / 8 + (intLen % 4 == 0 ? 0 : 1);

        write(intLen, c);
        write(byteLen, c);

        size_t b, i = 0;
        int mod4 = intLen % 4;
        for (b = 0; b < (mod4 == 0 ? byteLen : byteLen - 1); b++, i += 4) {
            c[b] = (data[i] << 6) | (data[i + 1] << 4) | (data[i + 2] << 2) | data[i + 3];
        }
        if (mod4 > 0) {
            if (mod4 == 1) {
                c[b] = (data[i] << 6);
            } else if (mod4 == 2) {
                c[b] = (data[i] << 6) | (data[i + 1] << 4);
            } else if (mod4 == 3) {
                c[b] = (data[i] << 6) | (data[i + 1] << 4) | (data[i + 2] << 2);
            }
        }
        c += byteLen;
    }


    std::vector<int> decode_int_2bits(const uchar *&c, size_t &remaining_length) {
        size_t byteLen, intLen;
        read(intLen, c, remaining_length);
        read(byteLen, c, remaining_length);
        std::vector<int> ints(intLen);
        size_t i = 0, b = 0;

        int mod4 = intLen % 4;
        for (; b < (mod4 == 0 ? byteLen : byteLen - 1); b++, i += 4) {
            ints[i] = (c[b] & 0xC0) >> 6;
            ints[i + 1] = (c[b] & 0x30) >> 4;
            ints[i + 2] = (c[b] & 0x0C) >> 2;
            ints[i + 3] = c[b] & 0x03;
        }
        if (mod4 > 0) {
            if (mod4 >= 1) {
                ints[i] = (c[b] & 0xC0) >> 6;
            }
            if (mod4 >= 2) {
                ints[i + 1] = (c[b] & 0x30) >> 4;
            }
            if (mod4 >= 3) {
                ints[i + 2] = (c[b] & 0x0C) >> 2;
            }
        }
        c += byteLen;
        remaining_length -= byteLen;
        return ints;
    }


    inline void encode_int_1bit(const std::vector<int> &data, uchar *&c) {

        size_t intLen = data.size();
        size_t byteLen = intLen / 8 + (intLen % 8 == 0 ? 0 : 1);

        write(intLen, c);
        write(byteLen, c);

        size_t b, i = 0;
        int mod8 = intLen % 8;
        for (b = 0; b < (mod8 == 0 ? byteLen : byteLen - 1); b++, i += 8) {
            c[b] = (data[i] << 7) | (data[i + 1] << 6) | (data[i + 2] << 5) | (data[i + 3] << 4)
                   | (data[i + 4] << 3) | (data[i + 5] << 2) | (data[i + 6] << 1) | (data[i + 7]);
        }
        if (mod8 > 0) {
            if (mod8 == 1) {
                c[b] = (data[i] << 7);
            } else if (mod8 == 2) {
                c[b] = (data[i] << 7) | (data[i + 1] << 6);
            } else if (mod8 == 3) {
                c[b] = (data[i] << 7) | (data[i + 1] << 6) | (data[i + 2] << 5);
            } else if (mod8 == 4) {
                c[b] = (data[i] << 7) | (data[i + 1] << 6) | (data[i + 2] << 5) | (data[i + 3] << 4);
            } else if (mod8 == 5) {
                c[b] = (data[i] << 7) | (data[i + 1] << 6) | (data[i + 2] << 5) | (data[i + 3] << 4)
                       | (data[i + 4] << 3);
            } else if (mod8 == 6) {
                c[b] = (data[i] << 7) | (data[i + 1] << 6) | (data[i + 2] << 5) | (data[i + 3] << 4)
                       | (data[i + 4] << 3) | (data[i + 5] << 2);
            } else if (mod8 == 7) {
                c[b] = (data[i] << 7) | (data[i + 1] << 6) | (data[i + 2] << 5) | (data[i + 3] << 4)
                       | (data[i + 4] << 3) | (data[i + 5] << 2) | (data[i + 6] << 1);
            }
        }
        c += byteLen;
    }

    std::vector<int> decode_int_1bit(const uchar *&c, size_t &remaining_length) {
        size_t byteLen, intLen;
        read(intLen, c, remaining_length);
        read(byteLen, c, remaining_length);
        std::vector<int> ints(intLen);
        size_t i = 0, b = 0;

        int mod8 = intLen % 8;
        for (; b < (mod8 == 0 ? byteLen : byteLen - 1); b++, i += 8) {
            ints[i] = (c[b] & 0x80) >> 7;
            ints[i + 1] = (c[b] & 0x40) >> 6;
            ints[i + 2] = (c[b] & 0x20) >> 5;
            ints[i + 3] = (c[b] & 0x10) >> 4;
            ints[i + 4] = (c[b] & 0x08) >> 3;
            ints[i + 5] = (c[b] & 0x04) >> 2;
            ints[i + 6] = (c[b] & 0x02) >> 1;
            ints[i + 7] = (c[b] & 0x01);
        }
        if (mod8 > 0) {
            if (mod8 >= 1) {
                ints[i] = (c[b] & 0x80) >> 7;
            }
            if (mod8 >= 2) {
                ints[i + 1] = (c[b] & 0x40) >> 6;
            }
            if (mod8 >= 3) {
                ints[i + 2] = (c[b] & 0x20) >> 5;
            }
            if (mod8 >= 4) {
                ints[i + 3] = (c[b] & 0x10) >> 4;
            }
            if (mod8 >= 5) {
                ints[i + 4] = (c[b] & 0x08) >> 3;
            }
            if (mod8 >= 6) {
                ints[i + 5] = (c[b] & 0x04) >> 2;
            }
            if (mod8 >= 7) {
                ints[i + 6] = (c[b] & 0x02) >> 1;
            }

        }
        c += byteLen;
        remaining_length -= byteLen;
        return ints;
    }

};
#endif //SZ3_BYTEUTIL_HPP