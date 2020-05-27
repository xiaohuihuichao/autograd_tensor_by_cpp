#include "Array.h"


bool isEqual(const float &_a, const float &_b)
{
    float sub = _a - _b;
    if (sub <= EPSINON && sub >= -EPSINON)
    {
        return true;
    }
    return false;
}

Array::Array(const int &_row, const int &_col, const float &_e) : m_row(_row), m_col(_col)
{
    if (_row <= 0 || _col <= 0)
    {
        throw Exception("Array::Array(const int &_row, const int &_col, const float &_e): the row and col of array must > 0.\n");
    }

    m_data = new float[m_row * m_col];
    for (int i = 0; i < m_row * m_col; ++i)
    {
        m_data[i] = _e;
    }
}

Array::Array(const Array &_arr) : m_row(_arr.m_row), m_col(_arr.m_col)
{
    if (_arr.m_data == NULL)
    {
        m_data = NULL;
        return;
    }

    m_data = new float[m_row * m_col];
    memcpy((void*)m_data, (void*)_arr.m_data, m_row * m_col * sizeof(float));
}

Array::~Array()
{
    release();
}

void Array::release()
{
    m_row = 0;
    m_col = 0;
    if (m_data != NULL)
    {
        delete m_data;
    }
    m_data = NULL;
}

Array diag(const int &_n, const float &_e)
{
    Array result(_n, _n);
    for (int i = 0; i < _n; ++i)
    {
        result.set(i, i, _e);
    }

    return result;
}

Array ones(const int &_row, const int &_col)
{
    return Array(_row, _col, 1);
}
Array ones_like(const Array &_arr)
{
    return Array(_arr.rows(), _arr.cols(), 1);
}
Array zeros(const int &_row, const int &_col)
{
    return Array(_row, _col, 0);
}
Array zeros_like(const Array &_arr)
{
    return Array(_arr.rows(), _arr.cols(), 0);
}

void Array::zero_()
{
    float *p = m_data;
    for (int i = 0; i < m_row * m_col; ++i, ++p)
    {
        *p = 0;
    }
}

float* Array::getRowPoint(const int &_row) const
{
    if (_row < 0 || _row >= m_row)
    {
        throw Exception("float* Array::getRowPoint(const int &_row) const: index is wrong.\n");
    }
    return m_data + _row * m_col;
}

float Array::get(const int &_row, const int &_col) const
{
    if (_row < 0 || _col < 0 || _row >= m_row || _col >= m_col)
    {
        throw Exception("float Array::get(const int &_row, const int &_col) const: index is wrong.\n");
    }

    float *pRow = getRowPoint(_row);
    return pRow[_col];
}

void Array::set(const int &_row, const int &_col, const float &_e)
{
    if (_row < 0 || _col < 0 || _row >= m_row || _col >= m_col)
    {
        throw Exception("void Array::set(const int &_row, const int &_col, const float &_e): index is wrong.\n");
    }

    float *pRow = getRowPoint(_row);
    pRow[_col] = _e;
}

void Array::swapRow(const int &_i, const int &_j)
{
    if (_i < 0 || _j < 0 || _i >= m_row || _j >= m_row)
    {
        throw Exception("void Array::swapRow(const int &_i, const int &_j): index is wrong.\n");
    }

    float tmp = 0;
    float *pI = getRowPoint(_i);
    float *pJ = getRowPoint(_j);
    for (int k = 0; k < m_col; ++k)
    {
        tmp = pI[k];
        pI[k] = pJ[k];
        pJ[k] = tmp;
    }
}

float Array::getDet() const
{
    if (m_col != m_row)
    {
        throw Exception("float Array::getDet() const: m_col != m_col.\n");
    }

    Array tmp = Array(*this);
    int ite = 0;
    float det = 1.;
    for (int i = 0; i < m_row; ++i)
    {
        if (isEqual(tmp(i, i), 0))
        {
            for (int j = i; j < m_row; ++j)
            {
                if (!isEqual(tmp(j, i), 0))
                {
                    tmp.swapRow(i, j);
                    ++ite;
                }
            }
        }

        for (int k = i + 1; k < m_row; ++k)
        {
            if (isEqual(tmp(i, i), 0))
            {
                return 0;
            }
            float scale = -1 * tmp(k, i) / tmp(i, i);
            float data = 0;
            for (int u = 0; u < m_row; ++u)
            {
                data = tmp(k, u) + tmp(i, u) * scale;
                tmp.set(k, u, data);
            }
        }
    }

    for (int i = 0; i < m_row; ++i)
    {
        det *= tmp(i, i);
    }

    if (1 == ite % 2)
    {
        det *= -1;
    }

    return det;
}

float Array::operator()(const int &_row, const int &_col) const
{
    return get(_row, _col);
}

Array& Array::operator=(const Array &_arr)
{
    if (_arr.m_data == m_data)
    {
        return *this;
    }

    if (NULL == m_data)
    {
        Array *r = new Array(_arr);
        this->m_row = r->m_row;
        this->m_col = r->m_col;
        this->m_data = r->m_data;
    }
    else if (m_row == _arr.m_row && m_col == _arr.m_col)
    {
        float *pDst = m_data;
        const float *pSrc = _arr.m_data;
        memcpy((void*)pDst, (void*)pSrc, m_row * m_col * sizeof(float));
    }
    else
    {
        throw Exception("Array& Array::operator=(const Array &_arr): size is not equal.\n");
    }

    return *this;
}

Array Array::inv_LU() const
{
    if (isEqual(getDet(), 0))
    {
        throw Exception("Array Array::inv_LU() const: determinant is equal to 0.\n");
    }

    int n = m_row;
    Array l(n, n), l_inverse(n, n);
    Array u(n, n), u_inverse(n, n);
    for (int i = 0; i < n; ++i)
    {
        l.set(i, i, 1);
        l_inverse.set(i, i, 1);
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            float s = 0;
            for (int k = 0; k < i; ++k)
            {
                s += l(i, k) * u(k, j);
            }
            u.set(i, j, get(i, j) - s);
        }
        for (int j = i + 1; j < n; ++j)
        {
            float s = 0;
            for (int k = 0; k < i; ++k)
            {
                s += l(j, k) * u(k, i);
            }
            float uii = u(i, i);
            if (isEqual(uii, 0))
            {
                throw Exception("Array Array::inv_LU() const: NAN.\n");
            }
            l.set(j, i, (get(j, i) - s) / uii);
        }
    }

    for (int i = 1; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            float s = 0;
            for (int k = 0; k < i; ++k)
            {
                s += l(i, k) * l_inverse(k, j);
            }
            l_inverse.set(i, j, -s);
        }
    }
    for (int i = 0; i < n; ++i)
    {
        float uii = u(i, i);
        if (isEqual(uii, 0))
        {
            throw Exception("Array Array::inv_LU() const: NAN.\n");
        }
        u_inverse.set(i, i, 1 / uii);
    }
    for (int i = 1; i < n; ++i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            float s = 0;
            for (int k = j + 1; k <= i; ++k)
            {
                s += u(j, k) * u_inverse(k, i);
            }
            float ujj = u(j, j);
            if (isEqual(ujj, 0))
            {
                throw Exception("Array Array::inv_LU() const: NAN.\n");
            }
            u_inverse.set(j, i, -s / ujj);
        }
    }

    return u_inverse * l_inverse;
}

Array Array::T() const
{
    if (m_row == 0 || m_col == 0)
    {
        return Array();
    }

    Array result(m_col, m_row);
    for (int r = 0; r < m_row; ++r)
    {
        float *pSrc = getRowPoint(r);
        for (int c = 0; c < m_col; ++c)
        {
            float *pDst = result.getRowPoint(c);
            pDst[r] = pSrc[c];
        }
    }

    return result;
}

Array Array::dot(const Array &_B) const
{
    return gemm_nn(*this, _B);
}


Array do_calc(const Array &_arr1, const Array &_arr2, calc_func func)
{
    Array result;
    if (_arr1.m_row == _arr2.m_row && _arr1.m_col == _arr2.m_col)
    {
        result = Array(_arr1);
        float *pDst = result.m_data;
        const float *pSrc2 = _arr2.m_data;
        for (int i = 0; i < _arr1.m_row * _arr1.m_col; ++i, ++pDst, ++pSrc2)
        {
            *pDst = (*func)(*pDst, *pSrc2);
        }
    }
    else if (_arr1.m_row == _arr2.m_row)
    {
        if (_arr1.m_col == 1)
        {
            result = Array(_arr2);
            float *pDst = result.m_data;
            const float *pSrc1 = _arr1.m_data;
            for (int i = 0; i < result.m_col; ++i)
            {
                for (int j = 0; j < result.m_col; ++j, ++pDst)
                {
                    *pDst = (*func)(pSrc1[i], *pDst);
                }
            }
        }
        else if (_arr2.m_col == 1)
        {
            result = Array(_arr1);
            float *pDst = result.m_data;
            const float *pSrc2 = _arr2.m_data;
            for (int i = 0; i < result.m_col; ++i)
            {
                for (int j = 0; j < result.m_col; ++j, ++pDst)
                {
                    *pDst = (*func)(pSrc2[i], *pDst);
                }
            }
        }
        else
        {
            throw Exception("Array do_calc(const Array &_arr1, const Array &_arr2, calc_func func): shape is not match.\n");
        }
    }
    else if (_arr1.m_col == _arr2.m_col)
    {
        if (_arr1.m_row == 1)
        {
            result = Array(_arr2);
            float *pDst = result.m_data;
            const float *pSrc1 = _arr1.m_data;
            for (int i = 0; i < result.m_col; ++i)
            {
                for (int j = 0; j < result.m_col; ++j, ++pDst)
                {
                    *pDst = (*func)(pSrc1[j], *pDst);
                }
            }
        }
        else if (_arr2.m_row == 1)
        {
            result = Array(_arr1);
            float *pDst = result.m_data;
            const float *pSrc2 = _arr2.m_data;
            for (int i = 0; i < result.m_col; ++i)
            {
                for (int j = 0; j < result.m_col; ++j, ++pDst)
                {
                    *pDst = (*func)(pSrc2[j], *pDst);
                }
            }
        }
        else
        {
            throw Exception("Array do_calc(const Array &_arr1, const Array &_arr2, calc_func func): shape is not match.\n");
        }
    }
    else if (_arr1.num_element() == 1)
    {
        return do_calc(*_arr1.m_data, _arr2, func);
    }
    else if (_arr2.num_element() == 1)
    {
        return do_calc(_arr1, *_arr2.m_data, func);
    }
    else
    {
        throw Exception("Array do_calc(const Array &_arr1, const Array &_arr2, calc_func func): shape is not match.\n");
    }

    return result;
}

Array do_calc(const Array &_arr, const float &_e, calc_func func)
{
    Array result(_arr);
    float *pDst = result.m_data;
    for (int i = 0; i < result.m_row * result.m_col; ++i, ++pDst)
    {
        *pDst = (*func)(*pDst, _e);
    }

    return result;
}
Array do_calc(const float &_e, const Array &_arr, calc_func func)
{
    Array result(_arr);
    float *pDst = result.m_data;
    for (int i = 0; i < result.m_row * result.m_col; ++i, ++pDst)
    {
        *pDst = (*func)(_e, *pDst);
    }

    return result;
}

Array& do_calc_(Array &_arr1, const Array &_arr2, calc_func func)
{
    if (_arr2.num_element() == 1)
    {
        return do_calc_(_arr1, *_arr2.m_data, func);
    }

    if (_arr1.m_col != _arr2.m_col)
    {
        throw Exception("Array& do_calc_(Array &_arr1, const Array &_arr2, calc_func func): col is not equal.\n");
    }

    float *pDst = _arr1.m_data;
    const float *pSrc = _arr2.m_data;

    if (_arr1.m_row == _arr2.m_row)
    {
        for (int i = 0; i < _arr1.m_row * _arr1.m_col; ++i, ++pDst, ++pSrc)
        {
            *pDst = (*func)(*pDst, *pSrc);
        }
    }
    else if (1 == _arr2.m_row)
    {
        for (int i = 0; i < _arr1.m_row; ++i)
        {
            for (int j = 0; j < _arr1.m_col; ++j, ++pDst)
            {
                *pDst = (*func)(*pDst, pSrc[j]);
            }
        }
    }
    
    return _arr1;
}
Array& do_calc_(Array &_arr, const float &_e, calc_func func)
{
    float *pDst = _arr.m_data;
    for (int i = 0; i < _arr.m_row * _arr.m_col; ++i, ++pDst)
    {
        *pDst = (*func)(*pDst, _e);
    }

    return _arr;
}

Array operator+(const Array &_arr1, const Array &_arr2)
{
    calc_func add_func = [=](float a, float b)-> float
    {
        return a + b;
    };

    return do_calc(_arr1, _arr2, add_func);
}
Array& operator+=(Array &_arr1, const Array &_arr2)
{
    calc_func add_func = [=](float a, float b)-> float
    {
        return a + b;
    };

    return do_calc_(_arr1, _arr2, add_func);
}
Array operator+(const Array &_arr, const float &_e)
{
    calc_func add_func = [=](float a, float b)-> float
    {
        return a + b;
    };

    return do_calc(_arr, _e, add_func);
}
Array& operator+=(Array &_arr, const float &_e)
{
    calc_func add_func = [=](float a, float b)-> float
    {
        return a + b;
    };

    return do_calc_(_arr, _e, add_func);
}
Array operator+(const float &_e, const Array &_arr)
{
    calc_func add_func = [=](float a, float b)-> float
    {
        return a + b;
    };

    return do_calc(_e, _arr, add_func);
}

Array operator*(const Array &_arr1, const Array &_arr2)
{
    calc_func mul_func = [=](float a, float b)-> float
    {
        return a * b;
    };

    return do_calc(_arr1, _arr2, mul_func);
}
Array& operator*=(Array &_arr1, const Array &_arr2)
{
    calc_func mul_func = [=](float a, float b)-> float
    {
        return a * b;
    };

    return do_calc_(_arr1, _arr2, mul_func);
}
Array operator*(const Array &_arr, const float &_e)
{
    calc_func mul_func = [=](float a, float b)-> float
    {
        return a * b;
    };

    return do_calc(_arr, _e, mul_func);
}
Array& operator*=(Array &_arr, const float &_e)
{
    calc_func mul_func = [=](float a, float b)-> float
    {
        return a * b;
    };

    return do_calc_(_arr, _e, mul_func);
}
Array operator*(const float &_e, const Array &_arr)
{
    calc_func mul_func = [=](float a, float b)-> float
    {
        return a * b;
    };

    return do_calc(_e, _arr, mul_func);
}

Array operator-(const Array &_arr)
{
    return -1 * _arr;
}
Array operator-(const Array &_arr1, const Array &_arr2)
{
    calc_func minus_func = [=](float a, float b)-> float
    {
        return a - b;
    };

    return do_calc(_arr1, _arr2, minus_func);
}
Array& operator-=(Array &_arr1, const Array &_arr2)
{
    calc_func minus_func = [=](float a, float b)-> float
    {
        return a - b;
    };

    return do_calc_(_arr1, _arr2, minus_func);
}
Array operator-(const Array &_arr, const float &_e)
{
    calc_func minus_func = [=](float a, float b)-> float
    {
        return a - b;
    };

    return do_calc(_arr, _e, minus_func);
}
Array& operator-=(Array &_arr, const float &_e)
{
    calc_func minus_func = [=](float a, float b)-> float
    {
        return a - b;
    };

    return do_calc_(_arr, _e, minus_func);
}
Array operator-(const float &_e, const Array &_arr)
{
    calc_func minus_func = [=](float a, float b)-> float
    {
        return a - b;
    };

    return do_calc(_e, _arr, minus_func);
}

Array operator/(const Array &_arr1, const Array &_arr2)
{
    calc_func div_func = [=](float a, float b)-> float
    {
        if (isEqual(b, 0))
        {
            throw Exception("Array operator/(const Array &_arr1, const Array &_arr2): / zero.\n");
        }
        return a / b;
    };

    return do_calc(_arr1, _arr2, div_func);
}
Array& operator/=(Array &_arr1, const Array &_arr2)
{
    calc_func div_func = [=](float a, float b)-> float
    {
        if (isEqual(b, 0))
        {
            throw Exception("Array& operator/=(Array &_arr1, const Array &_arr2): / zero.\n");
        }
        return a / b;
    };

    return do_calc_(_arr1, _arr2, div_func);
}
Array operator/(const Array &_arr, const float &_e)
{
    if (isEqual(_e, 0))
    {
        throw Exception("Array operator/(const Array &_arr, const float &_e): / zero.\n");
    }

    calc_func div_func = [=](float a, float b)-> float
    {
        return a / b;
    };
    
    return do_calc(_arr, _e, div_func);
}
Array& operator/=(Array &_arr, const float &_e)
{
    if (isEqual(_e, 0))
    {
        throw Exception("Array operator/(const Array &_arr, const float &_e): / zero.\n");
    }
    
    calc_func div_func = [=](float a, float b)-> float
    {
        return a / b;
    };
    
    return do_calc_(_arr, _e, div_func);
}
Array operator/(const float &_e, const Array &_arr)
{
    calc_func div_func = [=](float a, float b)-> float
    {
        if (isEqual(b, 0))
        {
            throw Exception("Array operator/(const Array &_arr, const float &_e): / zero.\n");
        }
        return a / b;
    };
    
    return do_calc(_e, _arr, div_func);
}

Array operator<(const Array &_arr1, const Array &_arr2)
{
    if (_arr1.m_row != _arr2.m_row || _arr1.m_col != _arr2.m_col)
    {
        throw Exception("Array operator<(const Array &_arr1, const Array &_arr2): shape not match.\n");
    }

    calc_func lt_func = [=](float a, float b)-> float
    {
        if (a < b)
        {
            return 1;
        }
        return 0;
    };

    return do_calc(_arr1, _arr2, lt_func);
}
Array operator<(const Array &_arr, const float &_e)
{
    calc_func lt_func = [=](float a, float b)-> float
    {
        if (a < b)
        {
            return 1;
        }
        return 0;
    };

    return do_calc(_arr, _e, lt_func);
}
Array operator<(const float &_e, const Array &_arr)
{
    calc_func lt_func = [=](float a, float b)-> float
    {
        if (a < b)
        {
            return 1;
        }
        return 0;
    };

    return do_calc(_e, _arr, lt_func);
}

Array operator>(const Array &_arr1, const Array &_arr2)
{
    if (_arr1.m_row != _arr2.m_row || _arr1.m_col != _arr2.m_col)
    {
        throw Exception("Array operator<(const Array &_arr1, const Array &_arr2): shape not match.\n");
    }

    calc_func lt_func = [=](float a, float b)-> float
    {
        if (a > b)
        {
            return 1;
        }
        return 0;
    };

    return do_calc(_arr1, _arr2, lt_func);
}
Array operator>(const Array &_arr, const float &_e)
{
    calc_func lt_func = [=](float a, float b)-> float
    {
        if (a > b)
        {
            return 1;
        }
        return 0;
    };

    return do_calc(_arr, _e, lt_func);
}
Array operator>(const float &_e, const Array &_arr)
{
    calc_func lt_func = [=](float a, float b)-> float
    {
        if (a > b)
        {
            return 1;
        }
        return 0;
    };

    return do_calc(_e, _arr, lt_func);
}

Array gemm_nn(const Array &_A, const Array &_B)
{
    if (_A.m_col != _B.m_row)
    {
        throw Exception("Array Array::gemm_nn(const Array &_A, const Array &_B): size is not match.\n");
    }

    Array result(_A.m_row, _B.m_col, 0);
    const float *pA = _A.m_data;
    float *pDst = result.m_data;
    for (int i = 0; i < _A.m_row; ++i, pDst += result.m_col)
    {
        for (int k = 0; k < _A.m_col; ++k, ++pA)
        {
            const float *pB = _B.getRowPoint(k);
            register float A_PART = *pA;
            for (int j = 0; j < _B.m_col; ++j)
            {
                pDst[j] += pB[j] * A_PART;
            }
        }
    }

    return result;
}
Array gemm_tn(const Array &_A, const Array &_B)
{
    if (_A.m_row != _B.m_row)
    {
        throw Exception("Array gemm_tn(const Array &_A, const Array &_B): size is not match.\n");
    }

    Array result(_A.m_col, _B.m_col, 0);
    float *pDst = result.m_data;
    for (int i = 0; i < _A.m_col; ++i, pDst += result.m_col)
    {
        for (int k = 0; k < _A.m_row; ++k)
        {
            register float A_PART = _A.get(k, i);
            const float *pB = _B.getRowPoint(k);
            for (int j = 0; j < _B.m_col; ++j)
            {
                pDst[j] += A_PART * pB[j];
            }
        }
    }

    return result;
}
Array gemm_nt(const Array &_A, const Array &_B)
{
    if (_A.m_col != _B.m_col)
    {
        throw Exception("Array gemm_nt(const Array &_A, const Array &_B): size is not match.\n");
    }

    Array result(_A.m_row, _B.m_row, 0);
    const float *pA = _A.m_data;
    float *pDst = result.m_data;
    for (int i = 0; i < _A.m_row; ++i, pA += _A.m_col)
    {
        for (int j = 0; j < _B.m_row; ++j, ++pDst)
        {
            register float sum = 0;
            const float *pB = _B.getRowPoint(j);
            for (int k = 0; k < _A.m_col; ++k)
            {
                sum += pA[k] * pB[k];
            }
            *pDst = sum;
        }
    }

    return result;
}
Array gemm_tt(const Array &_A, const Array &_B)
{
    if (_A.m_row != _B.m_col)
    {
        throw Exception("Array gemm_tt(const Array &_A, const Array &_B): size is not match.\n");
    }

    Array result(_A.m_col, _B.m_row, 0);
    float *pDst = result.m_data;
    for (int i = 0; i < _A.m_col; ++i)
    {
        for (int j = 0; j < _B.m_row; ++j, ++pDst)
        {
            register float sum = 0;
            const float *pB = _B.getRowPoint(j);
            for (int k = 0; k < _A.m_row; ++k)
            {
                sum += _A.get(k, i) * pB[k];
            }
            *pDst = sum;
        }
    }

    return result;
}
Array gemm(const Array &_A, const bool &_AT, const Array &_B, const bool &_BT)
{
    Array (*func) (const Array&, const Array&);

    if (!_AT && !_BT)
    {
        func = gemm_nn;
    }
    else if (_AT && !_BT)
    {
        func = gemm_tn;
    }
    else if (!_AT && _BT)
    {
        func = gemm_nt;
    }
    else
    {
        func = gemm_tt;
    }
    
    return (*func)(_A, _B);
}

Array& Array::element_do(element_func func)
{
    float *p = m_data;
    for (int i = 0; i < m_row * m_col; ++i, ++p)
    {
        *p = (*func)(*p);
    }

    return *this;
}

Array& Array::element_exp()
{
    return element_do(exp);
}
Array& Array::element_sqrt()
{
    return element_do(sqrt);
}
Array& Array::element_log2()
{
    return element_do(log);
}

double sigmoid(double x)
{
    return 1 / (1 + exp(x));
}
Array& Array::element_sigmoid()
{
    return element_do(sigmoid);
}

double tanh(double x)
{
    double exp_2x = exp(2*x);
    return (exp_2x - 1) / (exp_2x + 1);
}
Array& Array::element_tanh()
{
    return element_do(tanh);
}

double relu(double x)
{
    if (x < 0) return 0;

    return x;
}
Array& Array::element_relu()
{
    return element_do(relu);
}

Array& Array::element_pow(float n)
{
    float *p = m_data;
    for (int i = 0; i < m_row * m_col; ++i, ++p)
    {
        *p = pow(*p, n);
    }

    return *this;
}


Array Array::mean(int axis) const
{
    float num = 0;
    switch (axis)
    {
        case -1:
            num = m_row * m_col;
            break;
        case 0:
            num = m_row;
            break;
        case 1:
            num = m_col;
            break;
        default:
            throw Exception("axis is wrong.\n");
    }

    return sum(axis) / num;
}
Array Array::sum(int axis) const
{
    const float *p = m_data;
    Array result;
    float *pDst = NULL;

    switch (axis)
    {
        case -1:
            result = Array(1, 1, 0);
            pDst = result.m_data;
            for (int i = 0; i < m_row * m_col; ++i, ++p)
            {
                *pDst += *p;
            }
            break;
        case 0:
            result = Array(1, m_col, 0);
            pDst = result.m_data;
            for (int i = 0; i < m_row; ++i)
            {
                for (int j = 0; j < m_col; ++j, ++p)
                {
                    pDst[j] += *p;
                }
            }
            break;
        case 1:
            result = Array(m_row, 1, 0);
            pDst = result.m_data;
            for (int i = 0; i < m_row; ++i, ++pDst)
            {
                for (int j = 0; j < m_col; ++j, ++p)
                {
                    *pDst += *p;
                }
            }
            break;
        default:
            throw Exception("axis is wrong.\n");
    }

    return result;
}


std::ostream& operator<<(std::ostream &_out, const Array &_arr)
{
    float *p = _arr.m_data;
    _out << "Array ([\n";

    for (int i = 0; i < _arr.m_row; ++i)
    {
        _out << "  ";
        for (int j = 0; j < _arr.m_col; ++j, ++p)
        {
            _out << *p << ",\t";
        }
        _out << "\n";
    }
    _out << "], " << "shape: [" << _arr.m_row << ", " << _arr.m_col << "])\n";

    return _out;
}
