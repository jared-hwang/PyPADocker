
#ifndef bempp_cython_local_integrand_functors_hpp
#define bempp_cython_local_integrand_functors_hpp

#include "bempp/fiber/collection_of_3d_arrays.hpp"
#include "bempp/fiber/geometrical_data.hpp"
#include "bempp/fiber/simple_test_trial_integrand_functor.hpp"
#include "bempp/fiber/maxwell_3d_test_trial_integrand_functor.hpp"
#include "bempp/fiber/single_component_test_trial_integrand_functor.hpp"

namespace Fiber {



class LocalIntegrandFunctorBase
{
    public:

        virtual void addGeometricalDependencies(size_t &geomDeps) const = 0;

        inline virtual ~LocalIntegrandFunctorBase() {};

        virtual double evaluate(const ConstGeometricalDataSlice<double> &geomData,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& testValues,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& trialValues)
            const = 0;

};

template <typename Functor>
class ConcreteLocalIntegrandFunctor :
    public LocalIntegrandFunctorBase
{

    public:


        ConcreteLocalIntegrandFunctor(const Functor& functor) :
            m_functor(functor) {}

        virtual ~ConcreteLocalIntegrandFunctor() {};

        void addGeometricalDependencies(size_t &geomDeps) const override
        {
            m_functor.addGeometricalDependencies(geomDeps);
        }

        double evaluate(const ConstGeometricalDataSlice<double> &geomData,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& testValues,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& trialValues)
            const override {

                return m_functor.evaluate(geomData, testValues, trialValues);
            } 

    private:

        Functor m_functor;

};

class LocalIntegrandFunctorContainer
{
    public:
        typedef double CoordinateType;
        typedef double ResultType;
        typedef double BasisFunctionType;

        inline LocalIntegrandFunctorContainer(const shared_ptr<const LocalIntegrandFunctorBase>&
                functor) : m_functor(functor) {}

        inline void addGeometricalDependencies(size_t &geomDeps) const {

            m_functor->addGeometricalDependencies(geomDeps);

        }

        inline double evaluate(const ConstGeometricalDataSlice<double> &geomData,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& testValues,
                      const CollectionOf1dSlicesOfConst3dArrays<double>& trialValues)
            const {

                return m_functor->evaluate(geomData, testValues, trialValues);
            } 
        

    private:

        shared_ptr<const LocalIntegrandFunctorBase> m_functor;

};



inline LocalIntegrandFunctorContainer* simpleTestTrialIntegrandFunctor(){

    return new LocalIntegrandFunctorContainer(
                shared_ptr<LocalIntegrandFunctorBase>(
                    new ConcreteLocalIntegrandFunctor<SimpleTestTrialIntegrandFunctor<double, double>>(
                        SimpleTestTrialIntegrandFunctor<double, double>())
                    )
                );

}


inline LocalIntegrandFunctorContainer* maxwell3dTestTrialIntegrandFunctor(){

    return new LocalIntegrandFunctorContainer(
                shared_ptr<LocalIntegrandFunctorBase>(
                    new ConcreteLocalIntegrandFunctor<Maxwell3dTestTrialIntegrandFunctor<double, double>>(
                        Maxwell3dTestTrialIntegrandFunctor<double, double>())
                    )
                );

}

inline LocalIntegrandFunctorContainer* singleComponentTestTrialIntegrandFunctor(std::size_t testComponent, std::size_t trialComponent){

    return new LocalIntegrandFunctorContainer(
                shared_ptr<LocalIntegrandFunctorBase>(
                    new ConcreteLocalIntegrandFunctor<SingleComponentTestTrialIntegrandFunctor<double, double>>(
                        SingleComponentTestTrialIntegrandFunctor<double, double>(testComponent, trialComponent))
                    )
                );

}

}


#endif
