#include "Tensor.h"


Tensor::Tensor(const float &_e, const bool &_required_grad, const dependency &_dep): m_row(1), m_col(1), m_data(1, 1, _e), m_required_grad(_required_grad)
{
    if (m_required_grad)
    {
        m_grad = Array(1, 1, 0);
        dep = _dep;
    }
}
Tensor::Tensor(const Array &_arr, const bool &_required_grad, const dependency &_dep): m_row(_arr.rows()), m_col(_arr.cols()), m_data(_arr), m_required_grad(_required_grad)
{
    if (m_required_grad)
    {
        m_grad = Array(m_row, m_col, 0);
        dep = _dep;
    }
}
Tensor::Tensor(const int &_row, const int &_col, const float &_e, const bool &_required_grad): m_row(_row), m_col(_col), m_data(_row, _col, _e), m_required_grad(_required_grad) { if (m_required_grad) m_grad = Array(m_row, m_col, 0); }

Tensor::Tensor(const Tensor &_tensor, const bool &_required_grad): m_row(_tensor.m_row), m_col(_tensor.m_col), m_required_grad(_required_grad)
{
    m_data = Array(_tensor.m_data);
    if (m_required_grad) m_grad = Array(m_row, m_col, 0);
}

Tensor::~Tensor(){ release(); }

void Tensor::release()
{
    m_data.release();
    m_grad.release();
}

Tensor Tensor::dot(const Tensor &_t) const
{
    return tensor_dot(*this, _t);
}

Tensor Tensor::sum(int axis) const
{
    switch (axis)
    {
        case -1:
            return tensor_sum(*this);
    }
}

Tensor tensor_neg(const Tensor &_t)
{
    Array outData = -_t.m_data;
    bool required_grad = _t.m_required_grad;
    dependency dep;

    if (required_grad)
    {
        grad_fn func = [](const Array &grad)->Array
        {
            return -grad;
        };
        dep[&_t.m_data] = func;
    }

    return Tensor(outData, required_grad, dep);
}

Tensor tensor_dot(const Tensor &_t1, const Tensor &_t2)
{
    if (_t1.m_col != _t2.m_row)
    {
        throw Exception("Tensor tensor_dot(const Tensor &_t1, const Tensor &_t2): shape not match.\n");
    }

    Array outData = _t1.m_data.dot(_t2.m_data);
    bool required_grad = _t1.m_required_grad || _t2.m_required_grad;
    dependency dep;

    if (_t1.m_required_grad)
    {
        grad_fn func = [&_t2](const Array &grad)-> Array
        {
            return gemm(grad, false, _t2.m_data, true);
        };
        dep[&_t1.m_data] = func;
    }

    if (_t2.m_required_grad)
    {
        grad_fn func = [&_t1](const Array &grad)-> Array
        {
            return gemm(_t1.m_data, true, grad, false);
        };
        dep[&_t2.m_data] = func;
    }

    return Tensor(outData, required_grad, dep);
}

Tensor tensor_sum(const Tensor &_t)
{
    if (_t.num_element() == 0)
    {
        throw Exception("Tensor tensor_sum(const Tensor &_t): num of element is 0.\n");
    }

    Array outData = _t.m_data.sum();
    bool required_grad = _t.m_required_grad;
    dependency dep;

    if (required_grad)
    {
        grad_fn func = [&_t](const Array &grad)-> Array
        {
            return grad * ones_like(_t.m_data);
        };
        dep[&_t.m_data] = func;
    }

    return Tensor(outData, required_grad, dep);
}


int main()
{
    try
    {
        // Tensor bb(0, 0, 1, false);
        Tensor a(3, 4, 1, true), b(4, 5, 2, true);
        // Tensor c = tensor_dot(a, b);
        Tensor c = a.dot(b);
        dependency c_dep = c.get_dep();
        Array grad(3, 5, 1);
        for(dependency::iterator iter = c_dep.begin(); iter != c_dep.end(); iter++)
        {
            std::cout << iter->first << std::endl;
            std::cout << &a.data() << ", " << &b.data() << std::endl;
            std::cout << iter->second(grad) << std::endl;
        }

        Tensor sum_tensor(3, 8, 1, true);
        // Tensor s = tensor_sum(sum_tensor);
        Tensor s = sum_tensor.sum();
        c_dep = s.get_dep();
        Array grad_1(1, 1, 1);
        for(dependency::iterator iter = c_dep.begin(); iter != c_dep.end(); iter++)
        {
            std::cout << iter->first << std::endl;
            std::cout << &sum_tensor.data() << ", " << &s.data() << std::endl;
            std::cout << iter->second(grad_1) << std::endl;
        }
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
