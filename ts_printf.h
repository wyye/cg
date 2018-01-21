#ifndef TS_PRINTF_H
#define TS_PRINTF_H

#include <iostream>
#include <sstream>

namespace ts {

        inline void printf(std::ostream &stream, const char *format) {
                stream << format;
        }

        template<typename T, typename... Args>
        void printf(std::ostream &stream, const char *format, T val, Args... args) {
                const char *finish = format;
                while (*finish && *finish != '%') {
                        ++finish;
                }
                stream.write(format, finish - format);
                if (*finish == '%') {
                        stream << val;
                        printf(stream, finish + 1, args...);
                }
        }

        template<typename... Args>
        void printf(const char *format, Args... args) {
                printf(std::cout, format, args...);
        }

        template<typename... Args>
        std::string sprintf(const char *format, Args... args) {
                std::stringstream ss;
                printf(ss, format, args...);
                return ss.str();
        }

}

#endif
