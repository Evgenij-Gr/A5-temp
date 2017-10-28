#include "misc_defs.h"

int main(int argc, char* argv[])
{
    std::cout<<std::setprecision(8);
    // Test: finding events of harmonic oscillator
    {
        std::cout<<"[HARMONIC OSCILLATOR TESTS]"<<std::endl;
        bool testResult;
        int iterMax = 10;
        double atol = 1e-13;
        double rtol = 1e-13;
        double approximateReturnTime = 6.5;
        double crossingDirection = 1.0;
        int maxEvtsCapacity = iterMax +1;
        double evtTolerance = 1e-13;
        double timeSkip = 1e-3;
        HarmonicOscillator ho(1.0);
        HarmonicOscillatorEvent hoEvent;
        EventCalculator<HarmonicOscillator,
                HarmonicOscillatorEvent, HarmOscType> evcalc(atol, rtol,
                                                             approximateReturnTime,
                                                             crossingDirection,
                                                             iterMax,evtTolerance,
                                                             timeSkip);
        // ####################################################################
        {
            std::cout<<"@Starting on point not from a cross-section"<<std::endl;
            HarmOscType startPt = {3., 0.1};
            auto trajectoryStorage = evcalc.findEvents(ho, hoEvent, startPt);
            std::cout<<"Events: "<<std::endl;
            for (int i = 0; i <trajectoryStorage.size(); i++)
            {
                auto evt = trajectoryStorage[i];
                std::cout<<"Event[0] = ("<<evt[0]<<", "<<evt[1]<<")"<<std::endl;
            }
            {
                std::cout<<"[Test] Must have found at "
                           "least "<<(iterMax+1)<<" events"<<std::endl;
                testResult = (trajectoryStorage.size() >= maxEvtsCapacity );
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
            {
                std::cout<<"[Test] All X coordinates of "
                           "events should be negative"<<std::endl;
                testResult = std::all_of(trajectoryStorage.begin(),
                                         trajectoryStorage.end(),
                                         [](HarmOscType x){return (x[0]<0.0);});
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
            {
                std::cout<<"[Test] Starting point is not "
                           "the first point of trajectory"<<std::endl;
                testResult = (trajectoryStorage[0]!=startPt);
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
        }
        // ####################################################################
        {
            std::cout<<"@Starting on point from a cross-section"
                   " in wrong direction"<<std::endl;
            HarmOscType startPt = {3., 0.0};
            auto trajectoryStorage = evcalc.findEvents(ho, hoEvent, startPt);
            std::cout<<"Events: "<<std::endl;
            for (int i = 0; i <trajectoryStorage.size(); i++)
            {
                auto evt = trajectoryStorage[i];
                std::cout<<"Event[0] = ("<<evt[0]<<", "<<evt[1]<<")"<<std::endl;
            }
            {
                std::cout<<"[Test] Must have found at "
                       "least "<<(iterMax+1)<<" events"<<std::endl;
                testResult = (trajectoryStorage.size() >= maxEvtsCapacity );
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
            {
                std::cout<<"[Test] All X coordinates of "
                       "events should be negative"<<std::endl;
                testResult = std::all_of(trajectoryStorage.begin(),
                                     trajectoryStorage.end(),
                                     [](HarmOscType x){return (x[0]<0.0);});
                // there is a reason why this one doesn't fail in my problem
                // I put pt on a cross-section and choose crossing direction
                // depending on where it moves
                std::cout<<((testResult)?("PASSED"):("FAILED"))<<std::endl;
            }
        }
    }
    // Test: finding zeroes of sin-function
    // Test: finding zeroes of affine-function
    // Test: events count with direction
    // Test: finding fixed points
    // Test: eigenvalues precision
    return 0;
}
