#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <exception>
#include <string.h>
#include <math.h>


class Exception : public std::exception
{
private:
    const char *error_msg;

public:
    Exception(const char *_msg) { error_msg = _msg; }
    const char* what()
    {
        return error_msg;
    }
};


#define EPSINON 1e-6


typedef double (*element_func)(double);
typedef float (*calc_func)(float, float);

class Array
{
private:
    int m_row;
    int m_col;
    float *m_data;

public:
    Array() : m_row(0), m_col(0), m_data(NULL){}
    Array(const int &_row, const int &_col, const float &_e = 0);
    Array(const Array &_arr);
    ~Array();

    void release();

    int rows() const { return m_row; }
    int cols() const { return m_col; }
    int num_element() const { return m_row * m_col; }

    void zero_();

    float* getRowPoint(const int &_r) const;
    float get(const int &_row, const int &_col) const;
    void set(const int &_row, const int &_col, const float &_e);
    float getDet() const;
    Array inv_LU() const;

    Array T() const;

    float operator()(const int &_row, const int &_col) const;
    Array& operator=(const Array &_arr);

    Array dot(const Array &_B) const;

    friend Array do_calc(const Array &_arr1, const Array &_arr2, calc_func func);
    friend Array do_calc(const Array &_arr, const float &_e, calc_func func);
    friend Array do_calc(const float &_e, const Array &_arr, calc_func func);
    friend Array& do_calc_(Array &_arr1, const Array &_arr2, calc_func func);
    friend Array& do_calc_(Array &_arr, const float &_e, calc_func func);

    friend Array operator+(const Array &_arr1, const Array &_arr2);
    friend Array& operator+=(Array &_arr1, const Array &_arr2);
    friend Array operator+(const Array &_arr, const float &_e);
    friend Array& operator+=(Array &_arr, const float &_e);
    friend Array operator+(const float &_e, const Array &_arr);

    friend Array operator*(const Array &_arr1, const Array &_arr2);
    friend Array& operator*=(Array &_arr1, const Array &_arr2);
    friend Array operator*(const Array &_arr, const float &_e);
    friend Array& operator*=(Array &_arr, const float &_e);
    friend Array operator*(const float &_e, const Array &_arr);

    friend Array operator-(const Array &_arr);
    friend Array operator-(const Array &_arr1, const Array &_arr2);
    friend Array& operator-=(Array &_arr1, const Array &_arr2);
    friend Array operator-(const Array &_arr, const float &_e);
    friend Array& operator-=(Array &_arr, const float &_e);
    friend Array operator-(const float &_e, const Array &_arr);

    friend Array operator/(const Array &_arr1, const Array &_arr2);
    friend Array& operator/=(Array &_arr1, const Array &_arr2);
    friend Array operator/(const Array &_arr, const float &_e);
    friend Array& operator/=(Array &_arr, const float &_e);
    friend Array operator/(const float &_e, const Array &_arr);

    friend Array operator<(const Array &_arr1, const Array &_arr2);
    friend Array operator<(const Array &_arr, const float &_e);
    friend Array operator<(const float &_e, const Array &_arr);

    friend Array operator>(const Array &_arr1, const Array &_arr2);
    friend Array operator>(const Array &_arr, const float &_e);
    friend Array operator>(const float &_e, const Array &_arr);

    friend Array gemm_nn(const Array &_A, const Array &_B);
    friend Array gemm_tn(const Array &_A, const Array &_B);
    friend Array gemm_nt(const Array &_A, const Array &_B);
    friend Array gemm_tt(const Array &_A, const Array &_B);
    friend Array gemm(const Array &_A, const bool &_AT, const Array &_B, const bool &_BT);

    Array& element_exp();
    Array& element_sqrt();
    Array& element_log2();
    Array& element_sigmoid();
    Array& element_tanh();
    Array& element_relu();
    Array& element_pow(float n);

    Array sum(int axis=-1) const;
    Array mean(int axis=-1) const;

    friend std::ostream& operator<<(std::ostream &_out, const Array &_arr);
protected:
    void swapRow(const int &_i, const int &_j);
    Array& element_do(element_func func);
};

Array diag(const int &_n, const double &_e);
Array ones(const int &_row, const int &_col);
Array ones_like(const Array &_arr);
Array zeros(const int &_row, const int &_col);
Array zeros_like(const Array &_arr);


Array do_calc(const Array &_arr1, const Array &_arr2, calc_func func);
Array do_calc(const Array &_arr, const float &_e, calc_func func);
Array do_calc(const float &_e, const Array &_arr, calc_func func);
Array& do_calc_(Array &_arr1, const Array &_arr2, calc_func func);
Array& do_calc_(Array &_arr, const float &_e, calc_func func);

Array operator+(const Array &_arr1, const Array &_arr2);
Array& operator+=(Array &_arr1, const Array &_arr2);
Array operator+(const Array &_arr, const float &_e);
Array& operator+=(Array &_arr, const float &_e);
Array operator+(const float &_e, const Array &_arr);

Array operator*(const Array &_arr1, const Array &_arr2);
Array& operator*=(Array &_arr1, const Array &_arr2);
Array operator*(const Array &_arr, const float &_e);
Array& operator*=(Array &_arr, const float &_e);
Array operator*(const float &_e, const Array &_arr);

Array operator-(const Array &_arr);
Array operator-(const Array &_arr1, const Array &_arr2);
Array& operator-=(Array &_arr1, const Array &_arr2);
Array operator-(const Array &_arr, const float &_e);
Array& operator-=(Array &_arr, const float &_e);
Array operator-(const float &_e, const Array &_arr);

Array operator/(const Array &_arr1, const Array &_arr2);
Array& operator/=(Array &_arr1, const Array &_arr2);
Array operator/(const Array &_arr, const float &_e);
Array& operator/=(Array &_arr, const float &_e);
Array operator/(const float &_e, const Array &_arr);

Array operator<(const Array &_arr1, const Array &_arr2);
Array operator<(const Array &_arr, const float &_e);
Array operator<(const float &_e, const Array &_arr);

Array operator>(const Array &_arr1, const Array &_arr2);
Array operator>(const Array &_arr, const float &_e);
Array operator>(const float &_e, const Array &_arr);

Array gemm_nn(const Array &_A, const Array &_B);
Array gemm_tn(const Array &_A, const Array &_B);
Array gemm_nt(const Array &_A, const Array &_B);
Array gemm_tt(const Array &_A, const Array &_B);
Array gemm(const Array &_A, const bool &_AT, const Array &_B, const bool &_BT);
Array dot(const Array &_A, const Array &_B);

std::ostream& operator<<(std::ostream &_out, const Array &_arr);

#endif
