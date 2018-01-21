#include "fixedreal.h"

#include <iostream>

int numDigits(uint_fast64_t number)
{
    int result = 0;
    while (number) {
        number /= 10;
        result++;
    }
    return result;
}

int rtPow(int_fast64_t base, unsigned exp)
{
    int result = 1;
    while (exp) {
        exp -= 1;
        result *= base;
    }
    return result;
}

void test() {
    FixedReal<6> x = makeFixedReal(711, 4);
    std::cout << x << std::endl;

    x = makeFixedReal(723, 423);
    std::cout << x << std::endl;

    x = makeFixedReal(23);
    std::cout << x << std::endl;

    x = makeFixedReal(0, 423123);
    std::cout << x << std::endl;

    x = makeFixedReal(0, 0);
    std::cout << x << std::endl;

    FixedReal<6> y = makeFixedReal(2340.232234134);
    std::cout << y << std::endl;

    std::cout << y.toDouble() << std::endl;
    std::cout << y.toFloat() << std::endl;
    std::cout << "y = " << y.toString() << std::endl;

    x = makeFixedReal(21, 23423);
    std::cout << "x = " << x << std::endl;

    std::cout << "x + y = " << x + y << std::endl;
    std::cout << "-y = " << -y << std::endl;
    std::cout << "x + y = " << x + y << std::endl;
    std::cout << "x - y = " << x - y << std::endl;
    std::cout << "y / x = " << y / x << std::endl;
    std::cout << "y * x = " << y * x << std::endl;

    std::cout << makeFixedReal(-987654321987,654321) << std::endl;
    std::cout << makeFixedReal(907654321987,654821) - makeFixedReal(987654321987,654321) << std::endl;
}
