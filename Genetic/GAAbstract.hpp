#ifndef OptimT_GAABSTRACT_HPP
#define OptimT_GAABSTRACT_HPP

namespace OptimT {

template<class Args_t>
class parameterBody
{
public:
    parameterBody() {};
    virtual ~parameterBody() {};
    const Args_t & args() const {
        return _args;
    }

    void setArgs(const Args_t & a) {
        _args=a;
    }
protected:
Args_t _args;
};


template<>
class parameterBody<void>
{
public:
    parameterBody() {};
    virtual ~parameterBody() {};
};

}   //  OptimT

#endif  //  OptimT_GAABSTRACT_HPP