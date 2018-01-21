#ifndef LOGGER_H
#define LOGGER_H

#include <chrono>
#include <iomanip>
#include "ts_printf.h"

#ifndef LOG_LEVEL
#define LOG_LEVEL INFO
#endif

#define LOG_TYPES(F) \
        F(TRACE) \
        F(DEBUG) \
        F(INFO) \
        F(WARN) \
        F(ERROR)

#define DO_COMMA(OBJ) OBJ,
#define DO_ID(OBJ) OBJ

#define DO_DECL_FUNC(type) \
    template <typename... Args> \
    static void type(const char* format, Args... args) { \
                if (LogType::type >= LogType::LOG_LEVEL) { \
                        write(#type, format, args...); \
                } \
    }

class LOG {
private:
        enum class LogType {
                LOG_TYPES(DO_COMMA)
                OFF
        };
public:
        LOG_TYPES(DO_DECL_FUNC)

        class SECTION {
        public:
                template <typename... Args>
                SECTION(Args... args) :
                        m_message(ts::sprintf(args...))
                {
                        LOG::INFO("+++ %", m_message);
                }

                ~SECTION() {
                        LOG::INFO("=== %", m_message);
                }

        private:
                std::__cxx11::string m_message;
        };

private:
        template<typename... Args>
        static void write(const char *type, const char *format, Args... args) {
                std::chrono::time_point<std::chrono::_V2::system_clock> now;
                now = std::chrono::_V2::system_clock::now();
                auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
                time_t now_c = std::chrono::_V2::system_clock::to_time_t(now);
                std::ostringstream oss;
                oss << std::put_time(localtime(&now_c), "%F %T");
                oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
                std::__cxx11::string str = ts::sprintf(format, args...);
                ts::printf(std::cerr, "% [%] %\n", oss.str(), type[0], str);
        }
};

#undef LOG_TYPES
#undef DO_COMMA
#undef DO_ID
#undef DO_DECL_FUNC

#endif //LOGGER_LOGGER_H

