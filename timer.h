#ifndef TIMER_H
#define TIMER_H

#include <map>
#include <string>

class Timer {
  public:
    Timer() {
        started = false;
    }
    void X(std::string target) {
        double timeNow = glfwGetTime();
        if (!started) {
            started = true;
        } else {
            timeSpent[timingFor] += timeNow - previousTime;
        }
        previousTime = timeNow;
        timingFor = target;
    }
    ~Timer() {
        X("");
        for (std::map<std::string, double>::const_iterator it = timeSpent.begin();
             it != timeSpent.end(); ++it) {
            printf("%s: %f\n", it->first.c_str(), it->second);
        }
    }
  private:
    bool started;
    double previousTime;
    std::string timingFor;
    std::map<std::string, double> timeSpent;
};

#endif // TIMER_H
