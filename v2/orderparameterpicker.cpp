#include "orderparameterpicker.hpp"
#include "OrderParameter.hpp"
#include <stdexcept>

std::unique_ptr<OrderParameter> createOrderParameter(const std::string& name){
    //Insert your custom Order Parameter Reference here (Example to follow)
    // if (name =="OP"){
    //     return std::make_unique<MyOP>();
    // }

    throw std::runtime_error("Unknown CUSTOMOP '" + name +"'. Valid options: Biggest, AtomDistance");
}
