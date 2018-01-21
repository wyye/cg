#ifndef FIXEDREAL_H
#define FIXEDREAL_H

#include <cstdint>
#include <ostream>
#include <sstream>
#include <type_traits>
#include <iomanip>
#include <cmath>
#include <cstdlib>

template <int_fast64_t Base, unsigned Exp>
struct CTPow
{
    static const int_fast64_t value = Base * CTPow<Base, Exp - 1>::value;
};

template <int_fast64_t Base>
struct CTPow<Base, 0>
{
    static const int_fast64_t value = 1;
};

int rtPow(int_fast64_t base, unsigned exp);

int numDigits(uint_fast64_t number);

// 18 significant numbers (in base-10)
template <unsigned Frac>
struct FixedReal
{
    using Int = int_fast64_t;

    explicit FixedReal(Int internalValue) : value(internalValue) {}
    std::string toString() const {
        std::stringstream buffer;
        this->write(buffer);
        return buffer.str();
    }
    void write(std::ostream& os) const {
        if (sign()) {
            os << "-";
        }
        os << mag() << "." << std::setw(Frac) << std::setfill('0') << frac();
    }
    double toDouble() const { return double(value) / CTPow<10, Frac>::value; }
    float toFloat() const { return float(value) / CTPow<10, Frac>::value; }
    Int mag() const { return std::llabs(value / CTPow<10, Frac>::value); }
    Int frac() const { return std::llabs(value % CTPow<10, Frac>::value); }
    Int sign() const { return value < 0; }

    FixedReal<Frac> operator-() const { return FixedReal(-value); }

    template <Int Frac2>
    FixedReal<Frac2> changeFrac() const {
        if (Frac > Frac2) {
            return FixedReal<Frac2>(value / CTPow<10, Frac - Frac2>::value);
        } else {
            return FixedReal<Frac2>(value * CTPow<10, Frac2 - Frac>::value);
        }
    }

    Int getInternalValue() const { return value; }

private:
    Int value;
};

template <unsigned Frac = 6, typename T>
FixedReal<Frac> makeFixedReal(T mag) {
    return FixedReal<Frac>(mag * CTPow<10, Frac>::value);
}

template <unsigned Frac = 6, typename T>
FixedReal<Frac> makeFixedReal(T mag, uint_fast64_t frac) {
    static_assert(std::is_integral<T>::value, "magnitude should be integral");
    return FixedReal<Frac>(
                mag * CTPow<10, Frac>::value +
                frac * rtPow(10, Frac - numDigits(frac)));
}

template <unsigned Frac>
FixedReal<Frac> operator*(const FixedReal<Frac>& a, const FixedReal<Frac>& b) {
    return FixedReal<Frac>(a.getInternalValue() *
                           b.getInternalValue() /
                           CTPow<10, Frac>::value);
}

template <unsigned Frac>
FixedReal<Frac> operator/(const FixedReal<Frac>& a, const FixedReal<Frac>& b) {
    return FixedReal<Frac>(a.getInternalValue() *
                           CTPow<10, Frac>::value /
                           b.getInternalValue());
}

template <unsigned Frac>
FixedReal<Frac> operator+(const FixedReal<Frac>& a, const FixedReal<Frac>& b) {
    return FixedReal<Frac>(a.getInternalValue() + b.getInternalValue());
}

template <unsigned Frac>
FixedReal<Frac> operator-(const FixedReal<Frac>& a, const FixedReal<Frac>& b) {
    return FixedReal<Frac>(a.getInternalValue() - b.getInternalValue());
}

template <unsigned Frac>
std::ostream& operator<<(std::ostream& os, const FixedReal<Frac>& real) {
    real.write(os);
    return os;
}

void test();

#endif // FIXEDREAL_H
