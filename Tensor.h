#ifndef TENSOR_H
#define TENSOR_H

#include <map>
#include <functional>
#include "Array/Array.h"


typedef std::function<Array(const Array&)> grad_fn;
typedef std::map<const Array*, grad_fn> dependency;


class Tensor
{
private:
    int m_row;
    int m_col;
    Array m_data;
    
    bool m_required_grad;
    Array m_grad;
    dependency dep;

public:
    Tensor(): m_row(0), m_col(0), m_required_grad(false) {}
    Tensor(const float &_e, const bool &_required_grad=false, const dependency &_dep=dependency());
    Tensor(const Array &_arr, const bool &_required_grad=false, const dependency &_dep=dependency());
    Tensor(const int &_row, const int &_col, const float &_e=0, const bool &_required_grad=false);
    Tensor(const Tensor &_tensor, const bool &_required_grad=false);
    
    ~Tensor();
    
    void release();

    int rows() const { return m_row; }
    int cols() const { return m_col; }
    int num_element() const { return m_row * m_col; }
    const Array& data() const { return m_data; }


    ///********************************************
    dependency get_dep()
    {
        return dep;
    }
    //*********************************************/

    Tensor dot(const Tensor &_t) const;
    Tensor sum(int axis=-1) const;

    friend Tensor tensor_neg(const Tensor &_t);
    friend Tensor tensor_dot(const Tensor &_t1, const Tensor &_t2);
    friend Tensor tensor_sum(const Tensor &_t);
protected:
};

Tensor tensor_neg(const Tensor &_t);
Tensor tensor_dot(const Tensor &_t1, const Tensor &_t2);
Tensor tensor_sum(const Tensor &_t);

#endif
