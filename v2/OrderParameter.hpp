#pragma once
namespace LAMMPS_NS { class LAMMPS; } 
class OrderParameter{
    public:
    virtual ~OrderParameter() = default;
    virtual double update(LAMMPS_NS::LAMMPS* lmp) =0;
};