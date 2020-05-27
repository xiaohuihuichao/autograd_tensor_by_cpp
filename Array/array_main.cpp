#include <iostream>
#include "Array.h"

using namespace std;


void array_main()
{
    try
    {
        Array arr(1, 3, 2);
        Array cc_a(3, 3, 1);
        Array c_a(arr);
        Array ccc_a(cc_a);
        cc_a.set(1, 1, 100);
        cc_a.set(2, 2, 10);
        cc_a.set(0, 2, 7);

        // std::cout << cc_a.inv_LU() << std::endl;
        // std::cout << cc_a.getDet() << std::endl;

        // +
        std::cout << "cc_a: " << cc_a << "arr: " << arr << "cc_a + arr: " << cc_a + arr << std::endl;
        std::cout << "cc_a: " << cc_a << "cc_a + cc_a: " << cc_a + cc_a << std::endl;
        std::cout << "cc_a: " << cc_a << "e: " << 2 << "\ncc_a + e: " << cc_a + 2 << std::endl;
        std::cout << "e: " << 2 << "\ncc_a: " << cc_a << "e + cc_a: " << 2 + cc_a << std::endl;
        cc_a += arr;
        std::cout << "cc_a += arr: " << cc_a << std::endl;
        cc_a += cc_a;
        std::cout << "cc_a += cc_a: " << cc_a << std::endl;
        cc_a += 3;
        std::cout << "cc_a += 3: " << cc_a << std::endl;

        // -
        std::cout << "cc_a: " << cc_a << "arr: " << arr << "cc_a - arr: " << cc_a - arr << std::endl;
        std::cout << "cc_a: " << cc_a << "cc_a - cc_a: " << cc_a - cc_a << std::endl;
        std::cout << "cc_a: " << cc_a << "e: " << 2 << "\ncc_a - e: " << cc_a - 2 << std::endl;
        std::cout << "e: " << 2 << "\ncc_a: " << cc_a << "e - cc_a: " << 2 - cc_a << std::endl;
        cc_a -= arr;
        std::cout << "cc_a -= arr: " << cc_a << std::endl;
        cc_a -= cc_a;
        std::cout << "cc_a -= cc_a: " << cc_a << std::endl;
        cc_a -= 3;
        std::cout << "cc_a -= 3: " << cc_a << std::endl;

        // *
        std::cout << "cc_a: " << cc_a << "arr: " << arr << "cc_a * arr: " << cc_a * arr << std::endl;
        std::cout << "cc_a: " << cc_a << "cc_a * cc_a: " << cc_a * cc_a << std::endl;
        std::cout << "cc_a: " << cc_a << "e: " << 2 << "\ncc_a * e: " << cc_a * 2 << std::endl;
        std::cout << "e: " << 2 << "\ncc_a: " << cc_a << "e * cc_a: " << 2 * cc_a << std::endl;
        cc_a *= arr;
        std::cout << "cc_a *= arr: " << cc_a << std::endl;
        cc_a *= cc_a;
        std::cout << "cc_a *= cc_a: " << cc_a << std::endl;
        cc_a *= 3;
        std::cout << "cc_a *= 3: " << cc_a << std::endl;

        // /
        std::cout << "cc_a: " << cc_a << "arr: " << arr << "cc_a / arr: " << cc_a / arr << std::endl;
        std::cout << "cc_a: " << cc_a << "cc_a / cc_a: " << cc_a / cc_a << std::endl;
        std::cout << "cc_a: " << cc_a << "e: " << 2 << "\ncc_a / e: " << cc_a / 2 << std::endl;
        std::cout << "e: " << 2 << "\ncc_a: " << cc_a << "e / cc_a: " << 2 / cc_a << std::endl;
        cc_a /= arr;
        std::cout << "cc_a /= arr: " << cc_a << std::endl;
        cc_a /= cc_a;
        std::cout << "cc_a /= cc_a: " << cc_a << std::endl;
        cc_a /= 3;
        std::cout << "cc_a /= 3: " << cc_a << std::endl;

        // <
        Array arr_(3, 3, 1);
        Array out = cc_a < arr_;
        std::cout << "cc_a: " << cc_a << "arr_: " << arr_ << "cc_a < arr_: " << out << std::endl;
        out = cc_a < 2;
        std::cout << "cc_a: " << cc_a << "e: " << 2 << "\ncc_a < e: " << out << std::endl;
        out = 2 < cc_a;
        std::cout << "e: " << 2 << "\ncc_a: " << cc_a << "e < cc_a: " << out << std::endl;

        // >
        out = cc_a > arr_;
        std::cout << "cc_a: " << cc_a << "arr_: " << arr_ << "cc_a > arr_: " << out << std::endl;
        out = cc_a > 2;
        std::cout << "cc_a: " << cc_a << "e: " << 2 << "\ncc_a > e: " << out << std::endl;
        out = 2 > cc_a;
        std::cout << "e: " << 2 << "\ncc_a: " << cc_a << "e > cc_a: " << out << std::endl;

        Array A(3, 4);
        A.set(0, 0, 1);
        A.set(2, 1, 3);
        A.set(1, 3, 8);
        A.set(2, 2, 6);
        A.set(0, 3, 9);
        A.set(0, 2, 5);

        Array B(4, 2);
        B.set(3, 0, 1);
        B.set(2, 0, 3);
        B.set(1, 0, 8);
        B.set(2, 1, 6);
        B.set(0, 0, 9);
        B.set(0, 1, 5);
        std::cout << "A: " << A << "B: " << B << std::endl;
        std::cout << "gemm(A, B): " << gemm_nn(A, B) << std::endl;
        std::cout << "A.dot(B): " << A.dot(B) << std::endl;

        std::cout << "A.T(): " << A.T() << std::endl;
        std::cout << "gemm_tn(A.T(), B): " << gemm_tn(A.T(), B) << std::endl;
        std::cout << "gemm_nt(A, B.T()): " << gemm_nt(A, B.T()) << std::endl;
        std::cout << "gemm_tt(A.T(), B.T()): " << gemm_tt(A.T(), B.T()) << std::endl;
        
        std::cout << "gemm(A.T(), true, B.T(), true): " << gemm(A.T(), true, B.T(), true) << std::endl;
        std::cout << "gemm(A, false, B.T(), true): " << gemm(A, false, B.T(), true) << std::endl;
        std::cout << "gemm(A, false, B, false): " << gemm(A, false, B, false) << std::endl;
        std::cout << "gemm(A.T(), true, B, false): " << gemm(A.T(), true, B, false) << std::endl;

        std::cout << "A.sum: " << A.sum() << A.sum(0) << A.sum(1) << std::endl;
        std::cout << "A.mean: " << A.mean() << A.mean(0) << A.mean(1) << std::endl;
    }
    catch (Exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
    }
}


/*********************************
int main()
{
    array_main();
}
*********************************/
