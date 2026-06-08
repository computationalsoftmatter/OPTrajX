#pragma once
#include <memory>
#include <string>
#include "OrderParameter.hpp"

std::unique_ptr<OrderParameter> createOrderParameter(const std::string& name);