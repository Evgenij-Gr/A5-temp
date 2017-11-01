#ifndef DYN_UTILS
#define DYN_UTILS

#include <vector>
#include <cmath>

bool sameSign(double X, double Y, double zeroTol=1e-15)
{
    return ((X > zeroTol) && (Y > zeroTol)) || ((X < -zeroTol) && (Y < -zeroTol));
}

template<class Event, class Stepper, class System, class TStateType>
class EventObserver
{
private:
    System rhs;
    Event event;
    Stepper refiner;
public:
    std::vector< TStateType >& m_states;
private:
    double crossingDirection;
    double evtTolerance;
    bool isFirstCallToObserver;
    TStateType prevX;
    double prevT;
    double skipT;
    double prevEventTime;
public:
    EventObserver(System S, Event E, Stepper St, std::vector<TStateType> &states, double crDir, double evtTol, double skipTime) :
    rhs(S), event(E), refiner(St), m_states(states), crossingDirection(crDir), evtTolerance(evtTol), skipT(skipTime)
    {
//         prevEventTime = 0.0; // assume that we started exactly from cross-section? make it as a default value?
//         evtTolerance = 1e-15;
//        skipT = 1e-7;
        isFirstCallToObserver = true;
        // crossingDirection = 1.0;
    }
    void resetPreviousEventInformation()
    {
        isFirstCallToObserver = true;
    }
    void operator()(const TStateType &X, double t)
    {
//         std::cout<<"[integrator] X = ("<<X[0]<<", "<<X[1]<<", "<<X[2]<<", "<<X[3]<<") at t = "<<t<<std::endl;
        // this call establishes information about evts & times
        // in the start of integration
        if (isFirstCallToObserver)
        {
//             std::cout<<"[integrator] The first call to observer happened"<<std::endl;
            isFirstCallToObserver = false;
            if (std::abs(event(X)) < evtTolerance)
            {
//                 std::cout<<"[integrator] Started from an event point"<<std::endl;
//                 std::cout<<"[integrator] EVENT DETECTED"<<std::endl;
                m_states.push_back(X);
                prevEventTime = t; // basically it should be 0...
            }
            else
            {
//                 std::cout<<"[integrator] Started from a non-event point"<<std::endl;
                prevEventTime = t - skipT; // this is a trick, like -\infty for a first time when something occured if it hadn't
            }
            // prevX = X; // not needed really?
        }
        else
        {
            if ((t - prevEventTime) > skipT) // abs here? fix it? nah, assume that we always integrate in forward time
            {
                if (! sameSign(event(X), event(prevX)))
                {
//                     std::cout<<"[integrator] EVENT DETECTED"<<std::endl;
                    // event happened in between, let's determine when exactly
                    const int maxBisectionDepth = 70; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    int bisectionLevel = 0;
                    double dt = t - prevT;
                    refiner.initialize(prevX, prevT, dt);
                    auto timeRange = refiner.do_step(rhs);
                    double t0 = timeRange.first;
                    double t1 = timeRange.second;
                    double t_mid;
                    double event0 = event(prevX);
                    double event1 = event(X);
                    double dEvent = event1 - event0;
//                     std::cout<<"[integrator] event(X) = "<<event1<<std::endl;
//                     std::cout<<"[integrator] event(prevX) = "<<event0<<std::endl;
//                     std::cout<<"[integrator] dEvent = "<<dEvent<<std::endl;

                    TStateType x_mid;
                    // IN THIS CODE ALL CROSSINGS ARE FOUND NEED TO FIX IT
                    // SEEMS TO BE FIXED NOW
                    if ((crossingDirection * dEvent) > 1e-15)
                    {
                        while (( std::abs(dEvent) > std::abs(evtTolerance)) && (bisectionLevel < maxBisectionDepth))
                        {
                            ++bisectionLevel;
                            double t_mid = 0.5*(t0+t1);
                            // since we use dense output, do it other way
                            refiner.calc_state(t_mid, x_mid);
                            double event_mid = event(x_mid);
                            if (! sameSign(event0, event_mid))
                            {
                                t1 = t_mid;
                                event1 = event_mid;
                            }
                            else if (! sameSign(event_mid, event1))
                            {
                                t0 = t_mid;
                                event0 = event_mid;
                            }
                            dEvent = event1-event0;
                        }
                        // make the final step
                        t_mid = 0.5*(t0+t1);
                        refiner.calc_state(t_mid, x_mid);
                        // remember them
                        prevEventTime = t_mid;
                        // and write the result
//                         std::cout<<"[integrator] newEventTime = "<<t_mid<<std::endl;
//                         std::cout<<"[integrator] newEventX = ("<<x_mid[0]<<", "<<x_mid[1]<<", "<<x_mid[2]<<", "<<x_mid[3]<<")"<<std::endl;
                        event.adjustToSection(x_mid);
                        m_states.push_back( x_mid );
                    }
                }
            }
        }
        // this update happens always
        prevX = X;
        prevT = t;
    }
};

#endif
